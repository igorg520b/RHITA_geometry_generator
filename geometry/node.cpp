#include "node.h"


void icy::Node::Reset()
{
    x0.setZero();
    eqId = globId = -1;
    pinned = false;
    incident_faces.clear();
    surface = false;
}

void icy::Node::InitializeFromAnother(icy::Node *other)
{
    x0 = other->x0;
    pinned = other->pinned;
}


