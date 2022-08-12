#ifndef FL333_H
#define FL333_H

#include <vector>
#include <string>
#include <algorithm>

#include <gmsh.h>
#include <Eigen/Core>

#include <QString>
#include <QDir>
#include <QDebug>

#include "ObjectPool.h"
#include "node.h"

namespace icy { class Mesh; struct Node; struct Element; struct CohesiveZone;}

class icy::Mesh
{
public:
    Mesh();
    ~Mesh();
    Mesh& operator=(Mesh&) = delete;

    std::vector<icy::Node*> nodes;
    std::vector<icy::Element*> elems;
    std::vector<icy::CohesiveZone*> czs;

    icy::Node* AddNode();
    icy::Element* AddElement();
    icy::CohesiveZone* AddCZ();

    void Reset();

    void LoadMSH(const std::string &fileName, bool insertCZs);

    void ExportForAbaqus(std::string fileName, double czStrength, std::string jobName, std::string batchName,
                         double YoungsModulus, double czElasticity, double czEnergy,
                         bool rhitaSetup, double indenterRadius, double indenterDepth, double indentationRate,
                         double horizontalOffset, int nCPUs, int confinement);
    void RotateSample(double angleInDegrees);

private:
    constexpr static unsigned reserveConst = 100000;
    static SimplePool<Node> NodeFactory;
    static SimplePool<Element> ElementFactory;
    static SimplePool<CohesiveZone> CZFactory;
    void MarkIncidentFaces();

    // collision detection and contact
public:

    unsigned CountPinnedNodes() {return std::count_if(nodes.begin(),nodes.end(), [](Node* nd){return nd->pinned;});}

private:
    uint32_t FindClosestFace(Node *nd, Element *elem);
    static double dtn(Eigen::Vector3d pt, Eigen::Vector3d triangle[3]);
    void FaceCode2NodeCoords(uint32_t, Node *nds[3]) const;

};
#endif
