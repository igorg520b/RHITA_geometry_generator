#include "generator.h"
#include "gmsh.h"

#include "element2d.h"
#include "czinsertiontool2d.h"
#include "cohesivezone2d.h"

#include <spdlog/spdlog.h>
#include <fstream>
#include <iomanip>
#include <iostream>


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
//        if(nd->x0.y()==0 || nd->x0.x()==0 || abs(nd->x0.x()-blockLength)<1e-7) nd->group = 2;
        if(nd->x0.y()==0) nd->group = 2;
    }

    spdlog::info("blockLength {}; blockHeight {}", blockLength, blockHeight);
    spdlog::info("nds {}; elems {}; czs {}; grains {}", mesh2d.nodes.size(), mesh2d.elems.size(), mesh2d.czs.size(), dimTagsGrains.size());

//    indenterOffset = 0.05;
    notchOffset = 0;
    if(loadWithIndenter) CreatePyWithIndenter2D();

    spdlog::info("LoadFromFile done");
}



void Generator::CreatePyWithIndenter2D()
{
    spdlog::info("createVUEL {}", createVUEL);

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

    s << "elemType_bulk = mesh.ElemType(elemCode=CPS3, elemLibrary=STANDARD, secondOrderAccuracy=OFF,"
         "distortionControl=ON, lengthRatio=0.1, elemDeletion=ON)\n";

    bool hasCZs = mesh2d.czs.size()>0;

    if(hasCZs)
        s << "elemType_coh = mesh.ElemType(elemCode=COH2D4, elemLibrary=STANDARD,elemDeletion=ON)\n";


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

    // region - all nodes
    s << "region4allNodes = (";
    for(icy::Node2D *nd : mesh2d.nodes)
            s << "p.nodes["<<nd->globId<<":"<<nd->globId+1<<"],";
    s << ")\n";
    s << "p.Set(nodes=region4allNodes,name='Set4-all')\n";


    // create bulk material
    s << "mat1 = mdb.models['Model-1'].Material(name='Material-1-bulk')\n";
    s << "mat1.Density(table=((900.0, ), ))\n";
    s << "mat1.Elastic(table=((" << YoungsModulus << ", 0.3), ))\n";

    CreateCDP(s);

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
        s << "mat2.MaxsDamageInitiation(table=((" << czsStrength << "," << czsStrength*2 << "," << czsStrength*2 << "), ))\n";
        s << "mat2.maxsDamageInitiation.DamageEvolution(type=ENERGY, table=((" << czEnergy << ", ), ))\n";
        s << "mat2.Elastic(type=TRACTION, table=((" << czElasticity << "," << czElasticity << "," << czElasticity << "), ))\n";

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
    xOffset -= interactionRadius*1.5;
    spdlog::info("xOffset {}; indenterOffset {}; notchOffset {}",xOffset, indenterOffset, notchOffset);
    //horizontalOffset += -sqrt(pow(indenterRadius,2)-pow(indenterRadius-indenterDepth,2))-1e-7;
    double yOffset = blockHeight-d+R;

    double zOffset = 0;
    s << "a1.Instance(name='Part-2-1', part=p2, dependent=ON)\n";
    // rotate indenter
    s << "a1.rotate(instanceList=('Part-2-1', ), axisPoint=(0.0, 0.0, 0.0)," << "axisDirection=(0.0, 0.0, 1.0), angle=45.0)\n";
    s << "a1.translate(instanceList=('Part-2-1', ), vector=("<< xOffset << ", " << yOffset << ", " << zOffset << "))\n";

    // create step
    s << "mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=" << timeToRun << ", improvedDtMethod=ON)\n";

    // create field output request
    s << "mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(numIntervals=" << nFrames <<
         ",variables=('S', 'SVAVG', 'PE', 'PEVAVG', 'PEEQ', 'PEEQVAVG', 'LE', "
             "'U', 'V', 'A', 'RF', 'CSTRESS', 'DAMAGEC', 'DAMAGET', 'DAMAGESHR', 'EVF', "
             "'STATUS'))\n";

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
    s << "mdb.models['Model-1'].ContactProperty('IntProp-2czs')\n";

    s << "mdb.models['Model-1'].interactionProperties['IntProp-2czs'].TangentialBehavior("
         "formulation=FRICTIONLESS)\n";
    s << "mdb.models['Model-1'].interactionProperties['IntProp-2czs'].NormalBehavior("
         "pressureOverclosure=HARD, allowSeparation=ON, "
         "constraintEnforcementMethod=DEFAULT)\n";

//    s << "mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior("
//        "dependencies=0, directionality=ISOTROPIC, elasticSlipStiffness=None, "
//        "formulation=PENALTY, fraction=0.005, maximumElasticSlip=FRACTION, "
//        "pressureDependency=OFF, shearStressLimit=None, slipRateDependency=OFF, "
//        "table=((0.1, ), ), temperatureDependency=OFF)\n";

    s << "mdb.models['Model-1'].ContactProperty('IntProp-1')\n";
      s << "mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior("
        "formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, "
        "pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=(("
        "0.1, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, "
        "fraction=0.005, elasticSlipStiffness=None)\n";

       s << "mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior("
        "pressureOverclosure=EXPONENTIAL, table=((20000000.0, 0.0), (0.0, " << interactionRadius << ")), "
        "maxStiffness=None, constraintEnforcementMethod=DEFAULT)\n";


    // create interaction itself
//    s << "mdb.models['Model-1'].ContactExp(name='Int-1', createStepName='Step-1')\n";
//    s << "mdb.models['Model-1'].interactions['Int-1'].includedPairs.setValuesInStep("
//         "stepName='Step-1', useAllstar=ON)\n";
//    s << "mdb.models['Model-1'].interactions['Int-1'].contactPropertyAssignments.appendInStep("
//         "stepName='Step-1', assignments=((GLOBAL, SELF, 'IntProp-2czs'), ))\n";


    // additional contact between surface and nodes
    s << "mdb.models['Model-1'].SurfaceToSurfaceContactExp(clearanceRegion=None,"
        "createStepName='Initial', datumAxis=None, initialClearance=OMIT, "
        "interactionProperty='IntProp-1', main="
        "a1.surfaces['Surf-1'], "
        "mechanicalConstraint=KINEMATIC, name='Int-2', secondary="
        "a1.sets['MyPart1-1.Set4-all'], sliding=FINITE)\n";


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

void Generator::CreateCDP(std::ofstream &s)
{
    if(!createCDP) return;

    std::vector<double> CompressionStrainTable
    {
        0,
        0.000241232889562297,
        0.000354471872727266,
        0.000521980525252549,
        0.00071138618855219,
        0.000933655856565658,
        0.00118044919730642,
        0.00149661554545456,
        0.00194945956969697,
        0.00256248560673401,
        0.00333461398720537,
        0.00438538037710438,
        0.00896666666188552
    };

    std::vector<double> ConcreteCompressionHardening {
        27887324,
        37126760,
        38535212,
        39718308,
        40000000,
        39887324,
        39267604,
        38028168,
        36000000,
        32901408,
        28957748,
        24225352,
        400000
    };

    std::vector<double> ConcreteCompressionDamage {
        0,
        0,
        0,
        0,
        0,
        0.00281690000000001,
        0.0183099,
        0.0492958,
        0.1,
        0.1774648,
        0.2760563,
        0.3943662,
        0.99
    };



    // TENSILE PARAMETERS
    std::vector<double> &TensionStrainTable = CompressionStrainTable;

    std::vector<double> &ConcreteTensionStiffening = ConcreteCompressionHardening;

    std::vector<double> &ConcreteTensionDamage = ConcreteCompressionDamage;


    s << "mat1.ConcreteDamagedPlasticity(table=((40.0, 0.1, 1.16, 0.6667, 0.0), ))\n";

    s << std::setprecision(15);

    s << "mat1.concreteDamagedPlasticity.ConcreteCompressionHardening(table=(";
    for(int i=0;i<ConcreteCompressionHardening.size();i++)
    {
        const double stress = ConcreteCompressionHardening[i];
        const double strain = CompressionStrainTable[i];
        s << "(" << stress << ","  << strain << ")";
        if(i!=ConcreteCompressionHardening.size()-1) s << ",";
    }
    s << "))\n";

    s << "mat1.concreteDamagedPlasticity.ConcreteTensionStiffening(table=(";
    for(int i=0;i<ConcreteTensionStiffening.size();i++)
    {
        const double stress = ConcreteTensionStiffening[i];
        const double strain = TensionStrainTable[i];
        s << "(" << stress << ","  << strain << ")";
        if(i!=ConcreteTensionStiffening.size()-1) s << ",";
    }
    s << "))\n";

    s << "mat1.concreteDamagedPlasticity.ConcreteCompressionDamage(table=(";
    for(int i=0;i<ConcreteCompressionDamage.size();i++)
    {
        const double damage = ConcreteCompressionDamage[i];
        const double strain = CompressionStrainTable[i];
        s << "(" << damage << ","  << strain << ")";
        if(i!=ConcreteCompressionDamage.size()-1) s << ",";
    }
    s << "))\n";

    s << "mat1.concreteDamagedPlasticity.ConcreteTensionDamage(table=(";
    for(int i=0;i<ConcreteTensionDamage.size();i++)
    {
        const double damage = ConcreteTensionDamage[i];
        const double strain = TensionStrainTable[i];
        s << "(" << damage << ","  << strain << ")";
        if(i!=ConcreteTensionDamage.size()-1) s << ",";
    }
    s << "))\n";

    //*Concrete Failure,TYPE=Strain
    //0.0128889,  0.00348889, 0.9,  0.99

    std::cout << "*Concrete Failure,TYPE=Strain\n";
    std::cout << TensionStrainTable.back() << ",";
    std::cout << CompressionStrainTable.back() << ",";
    std::cout << ConcreteTensionDamage.back() << ",";
    std::cout << ConcreteCompressionDamage.back() << '\n';
}



void Generator::CreateTwoLayers(std::string MSHFileName)
{
    gmsh::initialize();

    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("indenter1");

    const double upperBlockHeight = indentationDepth*2;

    double upperBlockY = blockHeight - upperBlockHeight;

    int pt0 = gmsh::model::occ::addPoint(0, upperBlockY, 0);
    int pt1 = gmsh::model::occ::addPoint(0, upperBlockY+upperBlockHeight, 0);
    int pt2 = gmsh::model::occ::addPoint(blockLength, upperBlockY+upperBlockHeight, 0);
    int pt3 = gmsh::model::occ::addPoint(blockLength, upperBlockY, 0);

    int line1 = gmsh::model::occ::addLine(pt0,pt1);
    int line2 = gmsh::model::occ::addLine(pt1,pt2);
    int line3 = gmsh::model::occ::addLine(pt2,pt3);
    int line4 = gmsh::model::occ::addLine(pt3,pt0);

    int loopTag = gmsh::model::occ::addCurveLoop({line1,line2,line3,line4});
    gmsh::model::occ::addPlaneSurface({loopTag});

    gmsh::model::occ::synchronize();
    int groupTag2 = gmsh::model::addPhysicalGroup(1, {line4});

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

        icy::Node2D *nd = meshUpperBlock.AddNode();
        mtags[tag] = nd->globId;
        nd->x0 = Eigen::Vector2d(nodeCoords[i*3+0], nodeCoords[i*3+1]);
    }

    // mark node's groups
    nodeTags.clear();
    nodeCoords.clear();
    gmsh::model::mesh::getNodesForPhysicalGroup(1, groupTag2, nodeTags, nodeCoords);
    for(unsigned j=0;j<nodeTags.size();j++) meshUpperBlock.nodes[mtags[nodeTags[j]]]->group=2;
    spdlog::info("groupTag2 nodes {}",nodeTags.size());

    // get elements
    std::vector<std::size_t> trisTags, nodeTagsInTris;
    gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris);

    for(std::size_t i=0;i<trisTags.size();i++)
    {
        icy::Element2D *elem = meshUpperBlock.AddElement();
        elem->grainId = (int)i;
        for(int k=0;k<3;k++) elem->nds[k] = meshUpperBlock.nodes[mtags.at(nodeTagsInTris[i*3+k])];
    }


    for(icy::Element2D *elem : meshUpperBlock.elems) elem->Precompute();     // Dm matrix and volume

    spdlog::info("upper block: nodes {}; elems {}; czs {}", meshUpperBlock.nodes.size(), meshUpperBlock.elems.size(), meshUpperBlock.czs.size());
    meshUpperBlock.EvaluateMinMax();

    // UPPER BLOCK CREATED



    // CREATE LOWER BLOCK
    gmsh::clear();
    pt0 = gmsh::model::occ::addPoint(0, 0, 0);
    pt1 = gmsh::model::occ::addPoint(0, upperBlockY, 0);
    pt2 = gmsh::model::occ::addPoint(blockLength, upperBlockY, 0);
    pt3 = gmsh::model::occ::addPoint(blockLength, 0, 0);

    line1 = gmsh::model::occ::addLine(pt0,pt1);
    line2 = gmsh::model::occ::addLine(pt1,pt2);
    line3 = gmsh::model::occ::addLine(pt2,pt3);
    line4 = gmsh::model::occ::addLine(pt3,pt0);

    loopTag = gmsh::model::occ::addCurveLoop({line1,line2,line3,line4});
    gmsh::model::occ::addPlaneSurface({loopTag});

    gmsh::model::occ::synchronize();
    groupTag2 = gmsh::model::addPhysicalGroup(1, {line4});

    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMax", elemSize*3);
    gmsh::option::setNumber("Mesh.Algorithm", 5);
    gmsh::model::mesh::generate(2);

    // now load into MESH2D
    nodeTags.clear();
    nodeCoords.clear();
    parametricCoords.clear();
    mtags.clear(); // gmsh nodeTag -> sequential position in nodes[]

    // GET NODES
    gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);
    // set the size of the resulting nodes array
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in deformable");

        icy::Node2D *nd = meshLowerBlock.AddNode();
        mtags[tag] = nd->globId;
        nd->x0 = Eigen::Vector2d(nodeCoords[i*3+0], nodeCoords[i*3+1]);
    }

    // mark node's groups
    nodeTags.clear();
    nodeCoords.clear();
    gmsh::model::mesh::getNodesForPhysicalGroup(1, groupTag2, nodeTags, nodeCoords);
    for(unsigned j=0;j<nodeTags.size();j++) meshLowerBlock.nodes[mtags[nodeTags[j]]]->group=2;
    spdlog::info("groupTag2 nodes {}",nodeTags.size());

    // get elements
    trisTags.clear();
    nodeTagsInTris.clear();
    gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris);

    for(std::size_t i=0;i<trisTags.size();i++)
    {
        icy::Element2D *elem = meshLowerBlock.AddElement();
        elem->grainId = (int)i;
        for(int k=0;k<3;k++) elem->nds[k] = meshLowerBlock.nodes[mtags.at(nodeTagsInTris[i*3+k])];
    }

    // gmsh::write(outputFileName + ".msh");
    gmsh::finalize();

    for(icy::Element2D *elem : meshLowerBlock.elems) elem->Precompute();     // Dm matrix and volume

    spdlog::info("lower block: nodes {}; elems {}; czs {}", meshLowerBlock.nodes.size(), meshLowerBlock.elems.size(), meshLowerBlock.czs.size());
    meshLowerBlock.EvaluateMinMax();



    // CREATE .PY

    std::ofstream s;
    s.open(outputFileName+".py", std::ios_base::trunc|std::ios_base::out);
    s << std::setprecision(9);
    s << "from abaqus import *\n";
    s << "from abaqusConstants import *\n";
    s << "from caeModules import *\n";

    s << "import mesh\n";
    s << "import regionToolset\n";
    s << "import os\n";

    // types of elements
    s << "elemType_bulk = mesh.ElemType(elemCode=CPS3, elemLibrary=STANDARD, secondOrderAccuracy=OFF,"
         "distortionControl=ON, lengthRatio=0.1, elemDeletion=ON)\n";

    // part 1 (upper block)
    s << "p = mdb.models['Model-1'].Part(name='MyPart1', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)\n";

    for(icy::Node2D *nd : meshUpperBlock.nodes)
        s << "p.Node(coordinates=(" << nd->x0[0] << "," << nd->x0[1] << ",0))\n";

    s << "n = p.nodes\n";

    for(icy::Element2D *e : meshUpperBlock.elems)
        s << "p.Element(nodes=(n["<<e->nds[0]->globId<<"],n["<<e->nds[1]->globId<<
             "],n["<<e->nds[2]->globId<<"]), elemShape=TRI3)\n";

    // region1 - bulk elements
    s << "region1 = p.elements[0:" << meshUpperBlock.elems.size() << "]\n";
    s << "p.setElementType(regions=(region1,), elemTypes=(elemType_bulk,))\n";
    s << "p.Set(elements=(region1,), name='Set-1-elems')\n";


    // region - pinned nodes
    s << "region3pinned = (";
    for(icy::Node2D *nd : meshUpperBlock.nodes)
        if(nd->group==2)
            s << "p.nodes["<<nd->globId<<":"<<nd->globId+1<<"],";
    s << ")\n";
    s << "p.Set(nodes=region3pinned,name='Set3-pinned')\n";

    // region - all nodes
    s << "region4allNodes = (";
    for(icy::Node2D *nd : meshUpperBlock.nodes)
            s << "p.nodes["<<nd->globId<<":"<<nd->globId+1<<"],";
    s << ")\n";
    s << "p.Set(nodes=region4allNodes,name='Set4-all')\n";


    // bottom surface of the upper block
    std::vector<icy::Element2D*> f1elems, f2elems, f3elems;
    for(icy::Element2D *elem : meshUpperBlock.elems)
    {
        if(elem->nds[0]->x0.y() == meshUpperBlock.Ymin &&
                elem->nds[1]->x0.y() == meshUpperBlock.Ymin) f1elems.push_back(elem);
        if(elem->nds[1]->x0.y() == meshUpperBlock.Ymin &&
                elem->nds[2]->x0.y() == meshUpperBlock.Ymin) f2elems.push_back(elem);
        if(elem->nds[2]->x0.y() == meshUpperBlock.Ymin &&
                elem->nds[0]->x0.y() == meshUpperBlock.Ymin) f3elems.push_back(elem);
    }
    spdlog::info("f1elems {}; f2elems {}; f3elems {}", f1elems.size(), f2elems.size(), f3elems.size());

    s << "mdb.models['Model-1'].parts['MyPart1'].Surface(";
    if(!f1elems.empty())
    {
        s << "face1Elements=(";
        for(icy::Element2D *elem : f1elems) s << "p.elements[" << elem->elemId << ":" << elem->elemId+1 << "],";
        s << ")";
    }
    if(!f1elems.empty() && !f2elems.empty()) s << ",";
    if(!f2elems.empty())
    {
        s << "face2Elements=(";
        for(icy::Element2D *elem : f2elems) s << "p.elements[" << elem->elemId << ":" << elem->elemId+1 << "],";
        s << ")";
    }
    if(!(f2elems.empty() && f2elems.empty()) && !f3elems.empty()) s << ",";
    if(!f3elems.empty())
    {
        s << "face3Elements=(";
        for(icy::Element2D *elem : f3elems) s << "p.elements[" << elem->elemId << ":" << elem->elemId+1 << "],";
        s << ")";
    }
    s << ", name='Surf-1')\n";


    // part 2 (lower block)
    s << "p3lb = mdb.models['Model-1'].Part(name='PartLowerBlock', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)\n";

    for(icy::Node2D *nd : meshLowerBlock.nodes)
        s << "p3lb.Node(coordinates=(" << nd->x0[0] << "," << nd->x0[1] << ",0))\n";

    s << "n2 = p3lb.nodes\n";

    for(icy::Element2D *e : meshLowerBlock.elems)
        s << "p3lb.Element(nodes=(n2[" << e->nds[0]->globId << "],n2["<<e->nds[1]->globId <<
             "],n2["<<e->nds[2]->globId << "]), elemShape=TRI3)\n";

    // region1 - bulk elements
    s << "region11 = p3lb.elements[0:" << meshLowerBlock.elems.size() << "]\n";
    s << "p3lb.setElementType(regions=(region11,), elemTypes=(elemType_bulk,))\n";
    s << "p3lb.Set(elements=(region11,), name='Set-p2-elems')\n";


    // region - pinned nodes of the lower block
    s << "region4pinned = (";
    for(icy::Node2D *nd : meshLowerBlock.nodes)
        if(nd->group==2)
            s << "p3lb.nodes["<<nd->globId << ":"<<nd->globId+1 << "],";
    s << ")\n";
    s << "p3lb.Set(nodes=region4pinned,name='Set10-pinned-lower')\n";


    // top surface of the lower block
    f1elems.clear();
    f2elems.clear();
    f3elems.clear();
    for(icy::Element2D *elem : meshLowerBlock.elems)
    {
        if(elem->nds[0]->x0.y() == meshLowerBlock.Ymax &&
                elem->nds[1]->x0.y() == meshLowerBlock.Ymax) f1elems.push_back(elem);
        if(elem->nds[1]->x0.y() == meshLowerBlock.Ymax &&
                elem->nds[2]->x0.y() == meshLowerBlock.Ymax) f2elems.push_back(elem);
        if(elem->nds[2]->x0.y() == meshLowerBlock.Ymax &&
                elem->nds[0]->x0.y() == meshLowerBlock.Ymax) f3elems.push_back(elem);
    }
    spdlog::info("lower block: f1elems {}; f2elems {}; f3elems {}", f1elems.size(), f2elems.size(), f3elems.size());

    s << "p3lb.Surface(";
    if(!f1elems.empty())
    {
        s << "face1Elements=(";
        for(icy::Element2D *elem : f1elems) s << "p3lb.elements[" << elem->elemId << ":" << elem->elemId+1 << "],";
        s << ")";
    }
    if(!f1elems.empty() && !f2elems.empty()) s << ",";
    if(!f2elems.empty())
    {
        s << "face2Elements=(";
        for(icy::Element2D *elem : f2elems) s << "p3lb.elements[" << elem->elemId << ":" << elem->elemId+1 << "],";
        s << ")";
    }
    if(!(f2elems.empty() && f2elems.empty()) && !f3elems.empty()) s << ",";
    if(!f3elems.empty())
    {
        s << "face3Elements=(";
        for(icy::Element2D *elem : f3elems) s << "p3lb.elements[" << elem->elemId << ":" << elem->elemId+1 << "],";
        s << ")";
    }
    s << ", name='Surf-2')\n";


    // create bulk material
    s << "mat1 = mdb.models['Model-1'].Material(name='Material-1-bulk')\n";
    s << "mat1.Density(table=((900.0, ), ))\n";
    s << "mat1.Elastic(table=((" << YoungsModulus << ", 0.3), ))\n";

    CreateCDP(s);

    // bulk material without CDP
    s << "mat2 = mdb.models['Model-1'].Material(name='Material-2-bulk')\n";
    s << "mat2.Density(table=((900.0, ), ))\n";
    s << "mat2.Elastic(table=((" << YoungsModulus << ", 0.3), ))\n";

    // sections
    s << "mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1-CDP', "
         "material='Material-1-bulk', thickness=None)\n";
    s << "mdb.models['Model-1'].HomogeneousSolidSection(name='Section-2-bulk', "
         "material='Material-2-bulk', thickness=None)\n";

    // section assignments
    s << "region = p.sets['Set-1-elems']\n";
    s << "p.SectionAssignment(region=region, sectionName='Section-1-CDP', offset=0.0, "
         "offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)\n";

    // section - lower block
    s << "p2set = p3lb.sets['Set-p2-elems']\n";
    s << "p3lb.SectionAssignment(region=p2set, sectionName='Section-2-bulk', offset=0.0, "
         "offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)\n";


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
    s << "inst3 = a1.Instance(name='PartLowerBlock-1', part=p3lb, dependent=ON)\n";

    const double &R = indenterRadius;
    const double &d = indentationDepth;
    double b = sqrt(R*R - (R-d)*(R-d));
    double xOffset = -b;
    xOffset -= interactionRadius*1.3;
    spdlog::info("xOffset {}; indenterOffset {};",xOffset, indenterOffset);
    //horizontalOffset += -sqrt(pow(indenterRadius,2)-pow(indenterRadius-indenterDepth,2))-1e-7;
    double yOffset = blockHeight-d+R; //upperBlockHeight-d+R;

    double zOffset = 0;
    s << "a1.Instance(name='Part-2-1', part=p2, dependent=ON)\n";
    // rotate indenter
    s << "a1.rotate(instanceList=('Part-2-1', ), axisPoint=(0.0, 0.0, 0.0)," << "axisDirection=(0.0, 0.0, 1.0), angle=45.0)\n";
    s << "a1.translate(instanceList=('Part-2-1', ), vector=("<< xOffset << ", " << yOffset << ", " << zOffset << "))\n";

    // create step
    s << "mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=" << timeToRun << ", improvedDtMethod=ON)\n";

    // create field output request
    s << "mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(numIntervals=" << nFrames <<
         ",variables=('S', 'SVAVG', 'PE', 'PEVAVG', 'PEEQ', 'PEEQVAVG', 'LE', "
             "'U', 'V', 'A', 'RF', 'CSTRESS', 'DAMAGEC', 'DAMAGET', 'DAMAGESHR', 'EVF', "
             "'STATUS'))\n";

    // gravity load
    s << "mdb.models['Model-1'].Gravity(name='Load-1', createStepName='Step-1',comp2=-10.0, distributionType=UNIFORM, field='')\n";

    // BC - pinned nodes (lower block)
    s << "region = inst3.sets['Set10-pinned-lower']\n";
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
    s << "mdb.models['Model-1'].ContactProperty('IntProp-2czs')\n";

    s << "mdb.models['Model-1'].interactionProperties['IntProp-2czs'].TangentialBehavior("
         "formulation=FRICTIONLESS)\n";
    s << "mdb.models['Model-1'].interactionProperties['IntProp-2czs'].NormalBehavior("
         "pressureOverclosure=HARD, allowSeparation=ON, "
         "constraintEnforcementMethod=DEFAULT)\n";

    s << "mdb.models['Model-1'].ContactProperty('IntProp-1')\n";
      s << "mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior("
        "formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, "
        "pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=(("
        "0.1, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, "
        "fraction=0.005, elasticSlipStiffness=None)\n";

       s << "mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior("
        "pressureOverclosure=EXPONENTIAL, table=((20000000.0, 0.0), (0.0, " << interactionRadius << ")), "
        "maxStiffness=None, constraintEnforcementMethod=DEFAULT)\n";


    // contact between surface and nodes
    s << "mdb.models['Model-1'].SurfaceToSurfaceContactExp(clearanceRegion=None,"
        "createStepName='Initial', datumAxis=None, initialClearance=OMIT, "
        "interactionProperty='IntProp-1', main="
        "a1.surfaces['Surf-1'], "
        "mechanicalConstraint=KINEMATIC, name='Int-2', secondary="
        "a1.sets['MyPart1-1.Set4-all'], sliding=FINITE)\n";

    // tie constraint
    s << "mdb.models['Model-1'].Tie(adjust=ON, main=inst1.surfaces['Surf-1'],"
    "name='Constraint-2', positionToleranceMethod=COMPUTED,"
    "secondary=inst3.surfaces['Surf-2'], thickness=ON, tieRotations=ON)\n";


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
