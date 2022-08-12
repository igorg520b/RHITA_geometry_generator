#include "mesh.h"
#include "element.h"
#include "cohesivezone.h"
#include "czinsertiontool.h"

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


icy::SimplePool<icy::Node> icy::Mesh::NodeFactory(reserveConst);
icy::SimplePool<icy::Element> icy::Mesh::ElementFactory(reserveConst);
icy::SimplePool<icy::CohesiveZone> icy::Mesh::CZFactory(reserveConst);


icy::Mesh::Mesh()
{
}

icy::Mesh::~Mesh()
{
    Reset();
}

void icy::Mesh::Reset()
{
    NodeFactory.release(nodes);
    ElementFactory.release(elems);
    CZFactory.release(czs);
}



icy::Node* icy::Mesh::AddNode()
{
    Node* nd = NodeFactory.take();
    nd->Reset();
    nd->globId = (int)nodes.size();
    nodes.push_back(nd);
    return nd;
}

icy::Element* icy::Mesh::AddElement()
{
    Element* elem = ElementFactory.take();
    elem->Reset();
    elem->elemId = (int)elems.size();
    elems.push_back(elem);
    return elem;
}

icy::CohesiveZone* icy::Mesh::AddCZ()
{
    CohesiveZone *cz = CZFactory.take();
    cz->Reset();
    czs.push_back(cz);
    return cz;
}


void icy::Mesh::LoadMSH(const std::string &fileName, bool insertCZs)
{
    qDebug() << "LoadMSH";

    Reset();
    gmsh::clear();
    gmsh::open(fileName);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    std::unordered_map<std::size_t, std::size_t> mtags; // gmsh nodeTag -> sequential position in nodes[]

    // GET NODES
    gmsh::model::mesh::getNodesByElementType(4, nodeTags, nodeCoords, parametricCoords);


    NodeFactory.release(nodes);
    ElementFactory.release(elems);
    nodes.clear();
    elems.clear();

    // set the size of the resulting nodes array
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in deformable");

        Node *nd = AddNode();
        mtags[tag] = nd->globId;
        nd->x0 = Eigen::Vector3d(nodeCoords[i*3+0], nodeCoords[i*3+1], nodeCoords[i*3+2]);
    }

    // GET ELEMENTS - per grain (entity)
    std::vector<std::pair<int,int>> dimTagsGrains;
    gmsh::model::getEntities(dimTagsGrains,3);

    for(std::size_t j=0;j<dimTagsGrains.size();j++)
    {
        std::vector<std::size_t> tetraTags, nodeTagsInTetra;
        int entityTag = dimTagsGrains[j].second;
        gmsh::model::mesh::getElementsByType(4, tetraTags, nodeTagsInTetra,entityTag);

        for(std::size_t i=0;i<tetraTags.size();i++)
        {
            icy::Element *elem = AddElement();
            elem->grainId = (int)j;
            for(int k=0;k<4;k++) elem->nds[k] = nodes[mtags.at(nodeTagsInTetra[i*4+k])];
        }
    }

    std::vector<std::size_t> lineTags, nodeTagsInLines;
    gmsh::model::mesh::getElementsByType(1, lineTags, nodeTagsInLines);
    std::cout << "number of lines: " << lineTags.size() << std::endl;


    // center mesh
    for(Node *nd : nodes)
    {
        nd->pinned = nd->x0.z() < 1e-7;
    }


    CZInsertionTool czit;
    if(insertCZs) czit.InsertCZs(*this);

    for(icy::Element *elem : elems) elem->Precompute();     // Dm matrix and volume
    MarkIncidentFaces();

    gmsh::clear();
    qDebug() << "LoadMSH done";
}


void icy::Mesh::MarkIncidentFaces()
{
    // initialize incident faces information in nodes and elements

    // (1) find exposed faces
    std::map<std::tuple<int,int,int>,icy::CZInsertionTool::Facet> facets;
    for(icy::Element *e : elems)
    {
        for(int k=0;k<4;k++)
        {
            std::tuple<int,int,int> key = icy::CZInsertionTool::Facet::make_key(
                        e->nds[Element::fi[k][0]],
                    e->nds[Element::fi[k][1]],
                    e->nds[Element::fi[k][2]]);
            icy::CZInsertionTool::Facet facet;
            facet.key = key;
            facet.elems[0] = e;
            facet.facet_idx[0] = k;
            auto result = facets.insert({key,facet});
            if(result.second == false)
            {
                icy::CZInsertionTool::Facet &f = result.first->second;
                f.elems[1] = e;
                f.facet_idx[1] = k;
            }
        }
    }

    // (2) distribute exposed faces to their incident nodes
    for(Node *nd : nodes) { nd->incident_faces.clear(); nd->surface = false; }

    int count = 0;
    for(auto &kvp : facets)
    {
        icy::CZInsertionTool::Facet &f = kvp.second;
        if(f.elems[1] != nullptr) {
            continue;
        }

        Element *e = f.elems[0];
        unsigned k = f.facet_idx[0];

        uint32_t facet_code = (uint32_t)f.facet_idx[0] | (uint32_t)f.elems[0]->elemId << 2;
        for(int i=0;i<3;i++)
        {
            Node *nd = e->nds[Element::fi[k][i]];
            nd->incident_faces.push_back(facet_code);
            nd->surface = true;
        }
        count++;
    }

    // (3) per element - collect all exposed faces from element's 4 nodes
    for(Element *elem : elems)
    {
        auto &r = elem->elem_incident_faces;
        r.clear();
        for(int i=0;i<4;i++)
            for(uint32_t facet_code : elem->nds[i]->incident_faces)
                r.push_back(facet_code);
        std::sort(r.begin(),r.end());
        r.resize(std::distance(r.begin(),std::unique(r.begin(), r.end())));
    }
}











void icy::Mesh::RotateSample(double angleInDegrees)
{
    if(angleInDegrees == 0) return;
    double alpha = angleInDegrees*M_PI/180.;
    double cosA = std::cos(alpha);
    double sinA = std::sin(alpha);
    Eigen::Matrix3d R;
    R << cosA, -sinA, 0,
            sinA, cosA, 0,
            0, 0, 1;
    for(Node *nd : nodes)
    {
        nd->x0 = R*nd->x0;
    }
}
