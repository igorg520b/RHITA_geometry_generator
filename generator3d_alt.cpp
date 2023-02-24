#include "generator3d_alt.h"
#include "gmsh.h"

#include "element.h"
#include "czinsertiontool.h"
#include "cohesivezone.h"

#include <spdlog/spdlog.h>
#include <fstream>
#include <iomanip>
#include <iostream>




void Generator3DAlt::LoadFromFileWithCrop(std::string MSHFileName)
{
    spdlog::info("loading 3D from file {}",MSHFileName);

    gmsh::initialize();
    gmsh::open("msh3d\\" + MSHFileName);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    std::unordered_map<std::size_t, icy::Node*> mtags; // gmsh nodeTag -> node object

    gmsh::model::mesh::getNodesByElementType(4, nodeTags, nodeCoords, parametricCoords);

    // set the size of the resulting nodes array
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in deformable");
        icy::Node *nd = mesh.AddNode();
        mtags[tag] = nd;
        nd->x0 = Eigen::Vector3d(nodeCoords[i*3+0], nodeCoords[i*3+1], nodeCoords[i*3+2]);
    }


    // GET ELEMENTS - per grain (entity)
    std::vector<std::pair<int,int>> dimTagsGrains;
    gmsh::model::getEntities(dimTagsGrains,3);
    std::unordered_set<icy::Node*> used_nodes;

    for(std::size_t j=0;j<dimTagsGrains.size();j++)
    {
        std::vector<std::size_t> tetraTags, nodeTagsInTetra;
        int entityTag = dimTagsGrains[j].second;
        gmsh::model::mesh::getElementsByType(4, tetraTags, nodeTagsInTetra,entityTag);

        for(std::size_t i=0;i<tetraTags.size();i++)
        {
            icy::Node* nds[4];
            bool crop = false;
            for(int k=0;k<4;k++)
            {
                nds[k] = mtags.at(nodeTagsInTetra[i*4+k]);
                Eigen::Vector3d &vx = nds[k]->x0;
                if(makeCutout)
                {
                    double x = vx.x();
                    double y = vx.y();
                    if(isHalfSphere)
                    {
                        //                {&& x.z() < 0.151925 && x.x() > 5.65) crop = true;
                    }
                    else
                    {
                        // just a rectangular block
                        double xc = x-cutoutX;
                        double yc = y-(blockHeight+indenterRadius-indentationDepth);
                        if(xc*xc+yc*yc < indenterRadius*indenterRadius) crop = true;
                        if(y > blockHeight - indentationDepth && x < cutoutX) crop = true;
                    }
                }
            }
            if(crop) continue;

            icy::Element *elem = mesh.AddElement();
            elem->grainId = (int)j;
            elem->grainId = (int)j;
            for(int k=0;k<4;k++)
            {
                elem->nds[k] = nds[k];
                used_nodes.insert(nds[k]);
            }
        }
    }


    mesh.nodes.erase(std::remove_if(mesh.nodes.begin(), mesh.nodes.end(),
                                      [used_nodes](icy::Node *nd){ return used_nodes.count(nd)==0 ? true : false;}),
            mesh.nodes.end());
    for(int i=0;i<mesh.nodes.size();i++) mesh.nodes[i]->globId = i;

    gmsh::finalize();


    if(insertCZs)
    {
        icy::CZInsertionTool czit;
        if(insertCZs) czit.InsertCZs(mesh);
    }

    for(icy::Element *elem : mesh.elems) elem->Precompute();     // Dm matrix and volume


    // dimensions of the block
    auto it_x = std::max_element(mesh.nodes.begin(),mesh.nodes.end(),
                               [](icy::Node *nd1, icy::Node *nd2){return nd1->x0.x() < nd2->x0.x();});
    blockLength = (*it_x)->x0.x();

    auto it_y = std::max_element(mesh.nodes.begin(),mesh.nodes.end(),
                               [](icy::Node *nd1, icy::Node *nd2){return nd1->x0.y() < nd2->x0.y();});
    blockHeight = (*it_y)->x0.y();

    int count = 0;
    for(icy::Node *nd : mesh.nodes)
    {
        if(isHalfSphere)
        {}
        else
        {
            // recrangular block
            // mark attached layer
            if(nd->x0.y()==0) { nd->group = 2; continue; }
            // mark layer near the "indenter" for  force application
            double x = nd->x0.x();
            double y = nd->x0.y();
            double z = nd->x0.z();
            double xc = x-(cutoutX+deltaX);
            double yc = y-(blockHeight+indenterRadius-indentationDepth);
            if(xc*xc+yc*yc < indenterRadius*indenterRadius && z > 0.2 && z < 0.4)
            {
                nd->group = 55;
                count++;
            }
        }
    }
    const int nForcedNodes = std::count_if(mesh.nodes.begin(), mesh.nodes.end(),
                                           [](icy::Node* nd){return nd->group==55;});
    spdlog::info("forced nodes {}, {}",nForcedNodes,count);

    spdlog::info("blockLength {}; blockHeight {}", blockLength, blockHeight);
    spdlog::info("nds {}; elems {}; czs {}", mesh.nodes.size(), mesh.elems.size(), mesh.czs.size());

    CreatePy();
    spdlog::info("Generator3D::LoadFromFile done");
}



void Generator3DAlt::CreatePy()
{
    spdlog::info("Generator3D::CreatePy()");

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

    s << "p = mdb.models['Model-1'].Part(name='MyPart1', dimensionality=THREE_D, type=DEFORMABLE_BODY)\n";

    s << "print(\"importing " << mesh.nodes.size() << " nodes\")\n";

    for(int k=0;k<mesh.nodes.size();k++)
    {
        icy::Node* nd = mesh.nodes[k];
        if(k%(mesh.nodes.size()/100)==0)
            s << "print(\"" << 100*k/mesh.nodes.size() << "%\")\n";
        s << "p.Node(coordinates=(" << nd->x0[0] << "," << nd->x0[1] << "," << nd->x0[2] << "))\n";
    }

    s << "n = p.nodes\n";

    //print("nElems {0}; nCZS {1}; nTetra {2}".format(nElems,nCZS,nTetra))
    s << "print(\"importing " << mesh.elems.size() << " elements\")\n";
    for(int k=0;k<mesh.elems.size();k++)
    {
        icy::Element* e = mesh.elems[k];
        if(k%(mesh.elems.size()/100)==0)
            s << "print(\"" << 100*k/mesh.elems.size() << "%\")\n";
        s << "p.Element(nodes=(n["<<e->nds[0]->globId<<"],n["<<
             e->nds[1]->globId<<
             "],n["<<e->nds[3]->globId<<
             "],n["<<e->nds[2]->globId<<
             "]), elemShape=TET4)\n";
    }


    s << "print(\"importing " << mesh.czs.size() << " czs\")\n";
    for(int k=0;k<mesh.czs.size();k++)
    {
        icy::CohesiveZone* c = mesh.czs[k];
        if(k%(mesh.czs.size()/100)==0)
            s << "print(\"" << 100*k/mesh.czs.size() << "%\")\n";
        s << "p.Element(nodes=(n["<<c->nds[0]->globId<<
             "], n["<<c->nds[1]->globId<<
             "], n["<<c->nds[2]->globId<<
             "], n["<<c->nds[3]->globId<<
             "], n["<<c->nds[4]->globId<<
             "], n["<<c->nds[5]->globId<<"]), elemShape=WEDGE6)\n";
    }

    s << "elemType_bulk = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD, secondOrderAccuracy=OFF,"
         "distortionControl=ON, lengthRatio=0.1, elemDeletion=ON)\n";

    bool hasCZs = mesh.czs.size()>0;

    if(hasCZs)
    s << "elemType_coh = mesh.ElemType(elemCode=COH3D6, elemLibrary=STANDARD,elemDeletion=ON)\n";

    // region1 - bulk elements
    s << "region1 = p.elements[0:" << mesh.elems.size() << "]\n";
    s << "p.setElementType(regions=(region1,), elemTypes=(elemType_bulk,))\n";
    s << "p.Set(elements=(region1,), name='Set-1-elems')\n";

    if(hasCZs)
    {
        s << "region2cz = p.elements[" << mesh.elems.size() << ":" << mesh.elems.size() + mesh.czs.size() << "]\n";
        s << "p.setElementType(regions=(region2cz,), elemTypes=(elemType_coh,))\n";
        s << "p.Set(elements=(region2cz,), name='Set-2-czs')\n";
    }

    // region - pinned nodes

    s << "region3pinned = (";
    for(icy::Node *nd : mesh.nodes)
        if(nd->group==2)
            s << "p.nodes["<<nd->globId<<":"<<nd->globId+1<<"],";
    s << ")\n";
    s << "p.Set(nodes=region3pinned,name='Set3-pinned')\n";


    // region: "near" layer
    s << "regionForceLayer = (";
    for(icy::Node *nd : mesh.nodes)
    {
        if(nd->group==55)
            s << "p.nodes["<<nd->globId<<":"<<nd->globId+1<<"],";
    }
    s << ")\n";
    s << "p.Set(nodes=regionForceLayer,name='Set5-forceLayer')\n";


    // SURFACE
    //p = mdb.models['Model-1'].parts['MyPart1']
    //f = p.elements
    //face4Elements = f.getSequenceFromMask(mask=('[#0:447 #2000 ]', ), )
    //p.Surface(face4Elements=face4Elements, name='Surf-1')
    s << "f1Elems = (";
    for(int i=0;i<mesh.exteriorFacets[0].size();i++)
        s << "p.elements[" << mesh.exteriorFacets[0][i]->elemId << ":" << mesh.exteriorFacets[0][i]->elemId+1 << "],";
    s << ")\n";
    s << "f2Elems = (";
    for(int i=0;i<mesh.exteriorFacets[1].size();i++)
        s << "p.elements[" << mesh.exteriorFacets[1][i]->elemId << ":" << mesh.exteriorFacets[1][i]->elemId+1 << "],";
    s << ")\n";
    s << "f3Elems = (";
    for(int i=0;i<mesh.exteriorFacets[2].size();i++)
        s << "p.elements[" << mesh.exteriorFacets[2][i]->elemId << ":" << mesh.exteriorFacets[2][i]->elemId+1 << "],";
    s << ")\n";
    s << "f4Elems = (";
    for(int i=0;i<mesh.exteriorFacets[3].size();i++)
        s << "p.elements[" << mesh.exteriorFacets[3][i]->elemId << ":" << mesh.exteriorFacets[3][i]->elemId+1 << "],";
    s << ")\n";
    s << "p.Surface(face1Elements=f2Elems,face2Elements=f1Elems,face3Elements=f4Elems,face4Elements=f3Elems, name='Surf-1')\n";



    // create bulk material
    s << "mat1 = mdb.models['Model-1'].Material(name='Material-1-bulk')\n";
    s << "mat1.Density(table=((916.0, ), ))\n";
    s << "mat1.Elastic(table=((" << YoungsModulus << ", 0.3), ))\n";

    // cz material
    if(hasCZs)
    {
        s << "mat2 = mdb.models['Model-1'].Material(name='Material-2-czs')\n";
        s << "mat2.Density(table=((1.0, ), ))\n";
        s << "mat2.MaxsDamageInitiation(table=((" << czsStrength << "," << czsStrength*2 << "," << czsStrength*2 << "), ))\n";
        s << "mat2.maxsDamageInitiation.DamageEvolution(type=ENERGY, softening=EXPONENTIAL, table=((" << czEnergy << ", ), ))\n";
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

    // assembly
    s << "a1 = mdb.models['Model-1'].rootAssembly\n";
    s << "a1.DatumCsysByDefault(CARTESIAN)\n";

    // add and rotate main part
    s << "inst1 = a1.Instance(name='MyPart1-1', part=p, dependent=ON)\n";

    const double &R = indenterRadius;
    const double &d = indentationDepth;
    double b = sqrt(R*R - (R-d)*(R-d));
    double xOffset = -b;
    xOffset -= interactionRadius*1.5;
    spdlog::info("xOffset {}",xOffset);
    double yOffset = blockHeight-d+R;

    double zOffset = blockWidth/2;

    // create step
    s << "mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=" << timeToRun << ", improvedDtMethod=ON)\n";

    // create field output request
    s << "mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(numIntervals=" << nFrames <<
         ",variables=('S', 'SVAVG', 'PE', 'PEVAVG', 'PEEQ', 'PEEQVAVG', 'LE', "
             "'U', 'V', 'A', 'CSTRESS', 'DAMAGEC', 'DAMAGET', 'DAMAGESHR', 'EVF', "
             "'STATUS', 'SDEG'))\n";

    // gravity load
    s << "mdb.models['Model-1'].Gravity(name='Load-1', createStepName='Step-1',comp2=-10.0, distributionType=UNIFORM, field='')\n";

    // BC - pinned nodes
    s << "region = inst1.sets['Set3-pinned']\n";
    s << "mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Initial', region=region, localCsys=None)\n";


    // hydrostatic pressure
    s << "region = inst1.surfaces['Surf-1']\n";
    s << "mdb.models['Model-1'].SmoothStepAmplitude(name='Amp-2', timeSpan=STEP, data=((0.0, 0.0), (0.01, 1.0)))\n";

    s << "mdb.models['Model-1'].Pressure(name='Load-2', createStepName='Step-1', amplitude='Amp-2',"
        "region=region, distributionType=UNIFORM, field='', magnitude=1000000.0)\n";


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


