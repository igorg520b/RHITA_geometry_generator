#include "generator.h"
#include "gmsh.h"

#include "element2d.h"
#include "czinsertiontool2d.h"
#include "cohesivezone2d.h"

#include <spdlog/spdlog.h>
#include <fstream>
#include <iomanip>


void Generator::Generate()
{
    gmsh::initialize();

    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("indenter1");

    int pt0 = gmsh::model::occ::addPoint(0, blockHeight-indentationDepth,0);
    int pt1 = gmsh::model::occ::addPoint(0, blockHeight/2,0);
    int pt2 = gmsh::model::occ::addPoint(0,0,0);
    int pt3 = gmsh::model::occ::addPoint(blockLength,0,0);
    int pt4 = gmsh::model::occ::addPoint(blockLength, blockHeight/2, 0);
    int pt5 = gmsh::model::occ::addPoint(blockLength, blockHeight, 0);
    int pt6 = gmsh::model::occ::addPoint(notchOffset, blockHeight, 0);

    const double &R = indenterRadius;
    const double &d = indentationDepth;
    double b = sqrt(R*R - (R-d)*(R-d));

    int pt7 = gmsh::model::occ::addPoint(notchOffset-b,blockHeight-d,0);
    int ptC = gmsh::model::occ::addPoint(notchOffset-b, blockHeight-d+R,0);

    int arc =  gmsh::model::occ::addCircleArc(pt6,ptC,pt7);

    int line1 = gmsh::model::occ::addLine(pt7,pt0);
    int line2 = gmsh::model::occ::addLine(pt0,pt1);
    int line3 = gmsh::model::occ::addLine(pt1,pt2);
    int line4 = gmsh::model::occ::addLine(pt2,pt3);
    int line5 = gmsh::model::occ::addLine(pt3,pt4);
    int line6 = gmsh::model::occ::addLine(pt4,pt5);
    int line7 = gmsh::model::occ::addLine(pt5,pt6);

    int loopTag = gmsh::model::occ::addCurveLoop({arc,line1,line2,line3,line4,line5, line6, line7});
    gmsh::model::occ::addPlaneSurface({loopTag});

    gmsh::model::occ::synchronize();
    int groupTag1 = gmsh::model::addPhysicalGroup(1, {arc});
    int groupTag2 = gmsh::model::addPhysicalGroup(1, {line3, line4, line5});

    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMax", elemSize);
    gmsh::option::setNumber("Mesh.Algorithm", 5);
    gmsh::model::mesh::generate(2);

    // now load into MESH2D
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    std::unordered_map<std::size_t, std::size_t> mtags; // gmsh nodeTag -> sequential position in nodes[]

    // GET NODES
    gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);
    // set the size of the resulting nodes array
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in deformable");

        icy::Node2D *nd = mesh2d.AddNode();
        mtags[tag] = nd->globId;
        nd->x0 = Eigen::Vector2d(nodeCoords[i*3+0], nodeCoords[i*3+1]);
    }

    // mark node's groups
    nodeTags.clear();
    nodeCoords.clear();
    gmsh::model::mesh::getNodesForPhysicalGroup(1, groupTag1, nodeTags, nodeCoords);
    for(unsigned j=0;j<nodeTags.size();j++) mesh2d.nodes[mtags[nodeTags[j]]]->group=1;
    spdlog::info("groupTag1 nodes {}",nodeTags.size());

    nodeTags.clear();
    nodeCoords.clear();
    gmsh::model::mesh::getNodesForPhysicalGroup(1, groupTag2, nodeTags, nodeCoords);
    for(unsigned j=0;j<nodeTags.size();j++) mesh2d.nodes[mtags[nodeTags[j]]]->group=2;
    spdlog::info("groupTag2 nodes {}",nodeTags.size());

    // get elements
    std::vector<std::size_t> trisTags, nodeTagsInTris;
    gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris);

    for(std::size_t i=0;i<trisTags.size();i++)
    {
        icy::Element2D *elem = mesh2d.AddElement();
        elem->grainId = (int)i;
        for(int k=0;k<3;k++) elem->nds[k] = mesh2d.nodes[mtags.at(nodeTagsInTris[i*3+k])];
    }

    gmsh::write(outputFileName + ".msh");
    gmsh::finalize();



    // INSERT COHESIVE ZONES IF NEEDED

    icy::CZInsertionTool2D czit;
    if(insertCZs) czit.InsertCZs(mesh2d);

    for(icy::Element2D *elem : mesh2d.elems) elem->Precompute();     // Dm matrix and volume

    spdlog::info("nodes {}; elems {}; czs {}", mesh2d.nodes.size(), mesh2d.elems.size(), mesh2d.czs.size());

    if(loadWithIndenter) CreatePyWithIndenter2D();

}

void Generator::LoadFromFile(std::string MSHFileName)
{
    spdlog::info("loading from file {}",MSHFileName);

    gmsh::initialize();
    gmsh::open("msh\\" + MSHFileName);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    std::unordered_map<std::size_t, std::size_t> mtags; // gmsh nodeTag -> sequential position in nodes[]

    // GET NODES
    gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);

    // set the size of the resulting nodes array
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in deformable");
        icy::Node2D *nd = mesh2d.AddNode();
        mtags[tag] = nd->globId;
        nd->x0 = Eigen::Vector2d(nodeCoords[i*3+0], nodeCoords[i*3+1]);
    }

    // GET ELEMENTS - per grain (entity)
    std::vector<std::pair<int,int>> dimTagsGrains;
    gmsh::model::getEntities(dimTagsGrains,2);

    for(std::size_t j=0;j<dimTagsGrains.size();j++)
    {
        std::vector<std::size_t> trisTags, nodeTagsInTris;
        int entityTag = dimTagsGrains[j].second;
        gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris,entityTag);

        for(std::size_t i=0;i<trisTags.size();i++)
        {
            icy::Element2D *elem = mesh2d.AddElement();
            elem->grainId = (int)j;
            for(int k=0;k<3;k++) elem->nds[k] = mesh2d.nodes[mtags.at(nodeTagsInTris[i*3+k])];
        }
    }

    gmsh::finalize();

    if(insertCZs)
    {
        icy::CZInsertionTool2D czit;
        czit.InsertCZs(mesh2d);
    }

    for(icy::Element2D *elem : mesh2d.elems) elem->Precompute();     // Dm matrix and volume

    // dimensions of the block
    auto it_x = std::max_element(mesh2d.nodes.begin(),mesh2d.nodes.end(),
                               [](icy::Node2D *nd1, icy::Node2D *nd2){return nd1->x0.x() < nd2->x0.x();});
    blockLength = (*it_x)->x0.x();

    auto it_y = std::max_element(mesh2d.nodes.begin(),mesh2d.nodes.end(),
                               [](icy::Node2D *nd1, icy::Node2D *nd2){return nd1->x0.y() < nd2->x0.y();});
    blockHeight = (*it_y)->x0.y();

    for(icy::Node2D *nd : mesh2d.nodes)
    {
        if(nd->x0.y()==0 || nd->x0.x()==0 || abs(nd->x0.x()-blockLength)<1e-7) nd->group = 2;
    }

    spdlog::info("blockLength {}; blockHeight {}", blockLength, blockHeight);
    spdlog::info("nds {}; elems {}; czs {}; grains {}", mesh2d.nodes.size(), mesh2d.elems.size(), mesh2d.czs.size(), dimTagsGrains.size());

    indenterOffset = 0.05;
    if(loadWithIndenter) CreatePyWithIndenter2D();


    spdlog::info("LoadFromFile done");
}



void Generator::CreatePyWithIndenter2D()
{
    // SAVE as .py
    std::ofstream s;
    s.open(outputFileName+".py", std::ios_base::trunc|std::ios_base::out);
    s << std::setprecision(9);
    s << "from abaqus import *\n";
    s << "from abaqusConstants import *\n";
    s << "from caeModules import *\n";

    s << "import mesh\n";
    s << "import regionToolset\n";
    s << "import os\n";

    s << "p = mdb.models['Model-1'].Part(name='MyPart1', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)\n";

    for(icy::Node2D *nd : mesh2d.nodes)
        s << "p.Node(coordinates=(" << nd->x0[0] << "," << nd->x0[1] << ",0))\n";

    s << "n = p.nodes\n";

    for(icy::Element2D *e : mesh2d.elems)
        s << "p.Element(nodes=(n["<<e->nds[0]->globId<<"],n["<<e->nds[1]->globId<<
             "],n["<<e->nds[2]->globId<<"]), elemShape=TRI3)\n";

    for(icy::CohesiveZone2D *c : mesh2d.czs)
        s << "p.Element(nodes=(n["<<c->nds[0]->globId<<
             "], n["<<c->nds[1]->globId<<
             "], n["<<c->nds[3]->globId<<
             "], n["<<c->nds[2]->globId<<
             "]), elemShape=QUAD4)\n";

    s << "elemType_bulk = mesh.ElemType(elemCode=CPS3, elemLibrary=STANDARD, secondOrderAccuracy=OFF, distortionControl=DEFAULT)\n";

    bool hasCZs = mesh2d.czs.size()>0;
    if(hasCZs)
        s << "elemType_coh = mesh.ElemType(elemCode=COH2D4, elemLibrary=STANDARD)\n";


    // region1 - bulk elements
    s << "region1 = p.elements[0:" << mesh2d.elems.size() << "]\n";
    s << "p.setElementType(regions=(region1,), elemTypes=(elemType_bulk,))\n";
    s << "p.Set(elements=(region1,), name='Set-1-elems')\n";

    if(hasCZs)
    {
        s << "region2cz = p.elements[" << mesh2d.elems.size() << ":" << mesh2d.elems.size() + mesh2d.czs.size() << "]\n";
        s << "p.setElementType(regions=(region2cz,), elemTypes=(elemType_coh,))\n";
        s << "p.Set(elements=(region2cz,), name='Set-2-czs')\n";
    }


    // region - pinned nodes

    s << "region3pinned = (";
    for(icy::Node2D *nd : mesh2d.nodes)
        if(nd->group==2)
            s << "p.nodes["<<nd->globId<<":"<<nd->globId+1<<"],";
    s << ")\n";

    s << "p.Set(nodes=region3pinned,name='Set3-pinned')\n";


    // create bulk material
    s << "mat1 = mdb.models['Model-1'].Material(name='Material-1-bulk')\n";
    s << "mat1.Density(table=((900.0, ), ))\n";
    s << "mat1.Elastic(table=((" << YoungsModulus << ", 0.3), ))\n";

    if(plasticity)
    {
        s << "mat1.Plastic(scaleStress=None,table="
                "((50000.0, 0.0), (100000.0, 0.01), (1000000.0, 0.05), (10000000.0,0.1)))\n";
    }

    // cz material
    if(hasCZs)
    {
        s << "mat2 = mdb.models['Model-1'].Material(name='Material-2-czs')\n";
        s << "mat2.Density(table=((1.0, ), ))\n";
        s << "mat2.MaxsDamageInitiation(table=((" << czsStrength << "," << czsStrength/2 << "," << czsStrength/2 << "), ))\n";
        s << "mat2.maxsDamageInitiation.DamageEvolution(type=ENERGY, table=((" << czEnergy << ", ), ))\n";
        s << "mat2.Elastic(type=TRACTION, table=((" << czElasticity << "," << czElasticity/2 << "," << czElasticity/2 << "), ))\n";

        s << "mdb.models['Model-1'].CohesiveSection(name='Section-2-czs', "
             "material='Material-2-czs', response=TRACTION_SEPARATION, "
             "outOfPlaneThickness=None)\n";
    }

    // sections
    s << "mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1-bulk', "
         "material='Material-1-bulk', thickness=None)\n";

    // section assignments
    s << "region = p.sets['Set-1-elems']\n";
    s << "p.SectionAssignment(region=region, sectionName='Section-1-bulk', offset=0.0, "
         "offsetType=MIDDLE_SURFACE, offsetField='', "
         "thicknessAssignment=FROM_SECTION)\n";

    if(hasCZs)
    {
        s << "region = p.sets['Set-2-czs']\n";
        s << "p = mdb.models['Model-1'].parts['MyPart1']\n";
        s << "p.SectionAssignment(region=region, sectionName='Section-2-czs', offset=0.0, "
             "offsetType=MIDDLE_SURFACE, offsetField='', "
             "thicknessAssignment=FROM_SECTION)\n";
    }

    // indenter
    s << "s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2.0)\n";
    s << "g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n";
    s << "s.ArcByCenterEnds(center=(0.0, 0.0), point1=("<< -indenterRadius << ", 0.0), point2=(" << indenterRadius << ", -0.0125), direction=COUNTERCLOCKWISE)\n";

    s << "p2 = mdb.models['Model-1'].Part(name='Part-2', dimensionality=TWO_D_PLANAR, type=ANALYTIC_RIGID_SURFACE)\n";

    s << "p2.AnalyticRigidSurf2DPlanar(sketch=s)\n";


    s << "v1 = p2.vertices\n";

    s << "p2.ReferencePoint(point=p2.InterestingPoint(p2.edges[0], CENTER))\n";


    // assembly
    s << "a1 = mdb.models['Model-1'].rootAssembly\n";
    s << "a1.DatumCsysByDefault(CARTESIAN)\n";

    // add and rotate main part
    s << "inst1 = a1.Instance(name='MyPart1-1', part=p, dependent=ON)\n";

    const double &R = indenterRadius;
    const double &d = indentationDepth;
    double b = sqrt(R*R - (R-d)*(R-d));
    double xOffset = indenterOffset == 0 ? notchOffset-b : indenterOffset;
    double yOffset = blockHeight-d+R;
    double zOffset = 0;
    s << "a1.Instance(name='Part-2-1', part=p2, dependent=ON)\n";
    // rotate indenter
    s << "a1.rotate(instanceList=('Part-2-1', ), axisPoint=(0.0, 0.0, 0.0)," << "axisDirection=(0.0, 0.0, 1.0), angle=45.0)\n";
    s << "a1.translate(instanceList=('Part-2-1', ), vector=("<< xOffset << ", " << yOffset << ", " << zOffset << "))\n";

    // create step
    s << "mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=" << timeToRun << ", improvedDtMethod=ON)\n";

    // create field output request
    s << "mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(numIntervals=" << nFrames << ")\n";

    // gravity load
    s << "mdb.models['Model-1'].Gravity(name='Load-1', createStepName='Step-1',comp2=-10.0, distributionType=UNIFORM, field='')\n";

    // BC - pinned nodes
    s << "region = inst1.sets['Set3-pinned']\n";
    s << "mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Initial', region=region, localCsys=None)\n";

    // BC - moving indenter
    double xVelocity = indentationRate;
    double yVelocity = 0;
    s << "r1 = a1.instances['Part-2-1'].referencePoints\n";
    s << "refPoints1=(r1[2], )\n";
    s << "region = a1.Set(referencePoints=refPoints1, name='Set-1-indenterRP')\n";
    s << "mdb.models['Model-1'].VelocityBC(name='BC-2', createStepName='Step-1', "
         "region=region, v1="<< xVelocity << ", v2=" << yVelocity << ", v3=0.0, vr1=0.0, vr2=0.0, vr3=0.0, "
                                                                     "amplitude=UNSET, localCsys=None, distributionType=UNIFORM, fieldName='')\n";

    s << "ed1 = a1.instances['Part-2-1'].edges\n";
    s << "side2Edges1 = ed1[0:1]\n";
    s << "a1.Surface(name='Surf-1', side2Edges=side2Edges1)\n";

    s << "mdb.models['Model-1'].RigidBody(name='Constraint-1', refPointRegion=regionToolset.Region(referencePoints="
         "(a1.instances['Part-2-1'].referencePoints[2],)), surfaceRegion=a1.surfaces['Surf-1'])\n";


    // rigid body constraint

    // create interaction property
    s << "mdb.models['Model-1'].ContactProperty('IntProp-1')\n";
    s << "mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior("
         "formulation=FRICTIONLESS)\n";
    s << "mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior("
         "pressureOverclosure=HARD, allowSeparation=ON, "
         "constraintEnforcementMethod=DEFAULT)\n";

    // create interaction itself
    s << "mdb.models['Model-1'].ContactExp(name='Int-1', createStepName='Step-1')\n";
    s << "mdb.models['Model-1'].interactions['Int-1'].includedPairs.setValuesInStep("
         "stepName='Step-1', useAllstar=ON)\n";
    s << "mdb.models['Model-1'].interactions['Int-1'].contactPropertyAssignments.appendInStep("
         "stepName='Step-1', assignments=((GLOBAL, SELF, 'IntProp-1'), ))\n";


    // record indenter force

    s << "mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-1', name="
         "'H-Output-2', rebar=EXCLUDE, region="
         "mdb.models['Model-1'].rootAssembly.sets['Set-1-indenterRP'], sectionPoints="
         "DEFAULT, timeInterval=0.0001, variables=('RF1','RF2', ))\n";

    //create job
    s << "mdb.Job(name='" << outputFileName << "', model='Model-1', description='', type=ANALYSIS,"
                                        "atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,"
                                        "memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE,"
                                        "nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,"
                                        "contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='',"
                                        "resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains="
      <<numberOfCores<<","
                       "activateLoadBalancing=False, numThreadsPerMpiProcess=1,"
                       "multiprocessingMode=DEFAULT, numCpus="<<numberOfCores<<")\n";

    s.close();
}

