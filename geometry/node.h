#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing TBB headers
#ifndef NODE_H
#define NODE_H

#include <functional>
#include <vector>
#include <Eigen/Core>

namespace icy { struct Node; struct Element; struct CohesiveZone; }

struct icy::Node
{
    Node() { Reset(); incident_faces.reserve(8); }
    ~Node() = default;
    Node& operator=(Node&) = delete;

    void Reset();
    void InitializeFromAnother(icy::Node *other);

    int globId, eqId;               // id in fragment; id in mesh; id in freenode list;
    bool pinned;                    // the position of the node is set externally
    bool surface;

    Eigen::Vector3d x0;    // pos-initial, pos-current, velocity-current, pos-tentative
    Eigen::Vector3d F; // force applied via indenter

    std::vector<uint32_t> incident_faces; // list of faces that connect to this node; to highest bits indicate face # in element

};

#endif // NODE_H
#endif // Q_MOC_RUN
