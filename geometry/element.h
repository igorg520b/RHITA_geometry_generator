#ifndef ELEMENT123_H
#define ELEMENT123_H

#include "node.h"
#include <vector>
#include <tuple>
#include <bitset>
#include <algorithm>
#include <spdlog/spdlog.h>

namespace icy { struct Element; struct Node; }

struct icy::Element
{
    icy::Node* nds[4];          // initialized when the geometry is loaded or remeshed
    int grainId;
    int elemId;
    std::vector<uint32_t> elem_incident_faces; // faces that are incident to any of elem's nodes

    void Reset();
    void Precompute();

    double volume_initial;
    double strain_energy_density;   // (not multiplied by the volume)
    double principal_stress[3], max_shear_stress, hydrostatic_stress;

    Eigen::Matrix3d Dm, DmInv;  // reference shape matrix
    Eigen::Matrix3d CauchyStress, GreenStrain;

    void computeDm();


    int FacetToAbaqusFacetIndex(std::tuple<int,int,int>);



//    constexpr static int fi[4][3] {{1,2,3},{2,0,3},{3,0,1},{1,0,2}};
    // Abaqus face enumeration convention
    // Face 1: 1-2-3
    // Face 2: 1-4-2
    // Face 3: 2-4-3
    // Face 4: 3-4-1

    constexpr static int fi[4][3] {{0,1,2},{0,3,1},{1,3,2},{2,3,0}};

};

#endif // ELEMENT123_H
