#include "generator.h"
#include "gmsh.h"

#include "element2d.h"
#include "czinsertiontool2d.h"

#include <spdlog/spdlog.h>

Generator::Generator()
{

}

void Generator::Generate()
{
    gmsh::initialize();

    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("indenter1");

    int pt0 = gmsh::model::occ::addPoint(0, blockHeight-indentationDepth,0);
    int pt1 = gmsh::model::occ::addPoint(0,0,0);
    int pt2 = gmsh::model::occ::addPoint(blockLength,0,0);
    int pt3 = gmsh::model::occ::addPoint(blockLength, blockHeight, 0);
    int pt4 = gmsh::model::occ::addPoint(notchOffset, blockHeight, 0);
    spdlog::info("pt4 {},{}",notchOffset, blockHeight);

    const double &R = indenterRadius;
    const double &d = indentationDepth;
    double b = sqrt(R*R - (R-d)*(R-d));

    int pt5 = gmsh::model::occ::addPoint(notchOffset-b,blockHeight-d,0);
    int ptC = gmsh::model::occ::addPoint(notchOffset-b, blockHeight-d+R,0);

    int arc =  gmsh::model::occ::addCircleArc(pt4,ptC,pt5);

    int line1 = gmsh::model::occ::addLine(pt5,pt0);
    int line2 = gmsh::model::occ::addLine(pt0,pt1);
    int line3 = gmsh::model::occ::addLine(pt1,pt2);
    int line4 = gmsh::model::occ::addLine(pt2,pt3);
    int line5 = gmsh::model::occ::addLine(pt3,pt4);

    int loopTag = gmsh::model::occ::addCurveLoop({arc,line1,line2,line3,line4,line5});
    gmsh::model::occ::addPlaneSurface({loopTag});

    int groupTag1 = gmsh::model::addPhysicalGroup(1, {arc});
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMax", elemSize);
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

    // get elements
    std::vector<std::size_t> trisTags, nodeTagsInTris;
    gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris);

    for(std::size_t i=0;i<trisTags.size();i++)
    {
        icy::Element2D *elem = mesh2d.AddElement();
        elem->grainId = (int)i;
        for(int k=0;k<3;k++) elem->nds[k] = mesh2d.nodes[mtags.at(nodeTagsInTris[i*3+k])];
    }


    gmsh::write("test.msh");
    gmsh::finalize();

    icy::CZInsertionTool2D czit;
    if(insertCZs) czit.InsertCZs(mesh2d);

    for(icy::Element2D *elem : mesh2d.elems) elem->Precompute();     // Dm matrix and volume

    // SAVE as .py

}



