#include "mesh2d.h"
#include "element2d.h"
#include "cohesivezone2d.h"
#include "czinsertiontool2d.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <ios>
#include <iomanip>
#include <iostream>
#include <map>

#include <Eigen/Core>


icy::SimplePool<icy::Node2D> icy::Mesh2D::NodeFactory(reserveConst);
icy::SimplePool<icy::Element2D> icy::Mesh2D::ElementFactory(reserveConst);
icy::SimplePool<icy::CohesiveZone2D> icy::Mesh2D::CZFactory(reserveConst);


void icy::Mesh2D::Reset()
{
    NodeFactory.release(nodes);
    ElementFactory.release(elems);
    CZFactory.release(czs);
}

icy::Node2D* icy::Mesh2D::AddNode()
{
    Node2D* nd = NodeFactory.take();
    nd->Reset();
    nd->globId = (int)nodes.size();
    nodes.push_back(nd);
    return nd;
}

icy::Element2D* icy::Mesh2D::AddElement()
{
    Element2D* elem = ElementFactory.take();
    elem->Reset();
    elem->elemId = (int)elems.size();
    elems.push_back(elem);
    return elem;
}

icy::CohesiveZone2D* icy::Mesh2D::AddCZ()
{
    CohesiveZone2D *cz = CZFactory.take();
    cz->Reset();
    czs.push_back(cz);
    return cz;
}


void icy::Mesh2D::LoadMSH(const std::string &fileName, bool insertCZs)
{
    qDebug() << "LoadMSH - 2D";

    Reset();
    gmsh::clear();
    gmsh::open(fileName);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    std::unordered_map<std::size_t, std::size_t> mtags; // gmsh nodeTag -> sequential position in nodes[]

    // GET NODES
    gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);

    NodeFactory.release(nodes);
    ElementFactory.release(elems);
    CZFactory.release(czs);

    // set the size of the resulting nodes array
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in deformable");

        Node2D *nd = AddNode();
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
            icy::Element2D *elem = AddElement();
            elem->grainId = (int)j;
            for(int k=0;k<3;k++) elem->nds[k] = nodes[mtags.at(nodeTagsInTris[i*3+k])];
        }
    }



    CZInsertionTool2D czit;
    if(insertCZs) czit.InsertCZs(*this);

    for(icy::Element2D *elem : elems) elem->Precompute();     // Dm matrix and volume

    gmsh::clear();
    qDebug() << "LoadMSH 2D done";
}

