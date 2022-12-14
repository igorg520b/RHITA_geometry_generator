#ifndef ELEMENT123_H
#define ELEMENT123_H

#include "node.h"

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


    constexpr static int fi[4][3] {{1,2,3},{2,0,3},{3,0,1},{1,0,2}};

private:
    static Eigen::Matrix3d DDs[12]; // derivatives of Ds with respect to x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
    static Eigen::Matrix<double,12,12> consistentMassMatrix;

};

#endif // ELEMENT123_H
