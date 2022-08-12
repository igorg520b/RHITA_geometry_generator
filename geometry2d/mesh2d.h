#ifndef MESH2D_H
#define MESH2D_H

#include <vector>
#include <string>
#include <algorithm>

#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>

#include <gmsh.h>
#include <Eigen/Core>

#include <QString>
#include <QDir>
#include <QDebug>

#include "ConcurrentPool.h"
#include "node2d.h"

namespace icy { class Mesh2D; struct Node2D; struct Element2D; struct CohesiveZone2D;}

class icy::Mesh2D
{
public:
    Mesh2D() = default;
    ~Mesh2D() {Reset();};
    Mesh2D& operator=(Mesh2D&) = delete;

    std::vector<icy::Node2D*> nodes;
    std::vector<icy::Element2D*> elems;
    std::vector<icy::CohesiveZone2D*> czs;

    icy::Node2D* AddNode();
    icy::Element2D* AddElement();
    icy::CohesiveZone2D* AddCZ();

    void Reset();

    void LoadMSH(const std::string &fileName, bool insertCZs);

private:
    constexpr static unsigned reserveConst = 100000;
    static ConcurrentPool<Node2D> NodeFactory;
    static ConcurrentPool<Element2D> ElementFactory;
    static ConcurrentPool<CohesiveZone2D> CZFactory;
//    void MarkIncidentFaces();

};
#endif
