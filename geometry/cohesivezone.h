#ifndef COHESIVEZONE_H
#define COHESIVEZONE_H

#include "node.h"
#include <cstdint>

namespace icy {struct CohesiveZone; struct Node; struct Element;}

struct icy::CohesiveZone
{
    Element *elems[2];                  // each CZ connects two elements
    uint8_t faceIds[2];
    Node* nds[6];

    constexpr static int nQPts = 3;     // number of quadrature points (not the number of nodes!)
    double pmax[nQPts];                 // max normal separation reached
    double tmax[nQPts];                 // max tangential separation reached
    bool isActive;                      // if not active, CZ has failed
    bool isDamaged;

    void Reset();

private:
    bool tentative_contact, tentative_failed, tentative_damaged;
    double tentative_pmax_final, tentative_tmax_final;
    double tentative_pmax[nQPts], tentative_tmax[nQPts];

    constexpr static double epsilon = -1e-9;
    constexpr static double epsilon_abs = 1e-9;
    constexpr static double epsilon_fail_traction = 0.05; // if traction is <5% of max, CZ fails

    const static Eigen::Vector3d quadraturePoints[nQPts];
    const static Eigen::Matrix<double,18,3> B[nQPts];


};

#endif // COHESIVEZONE_H
