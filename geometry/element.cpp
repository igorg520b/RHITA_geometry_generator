#include "element.h"
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include <complex>
#include <cmath>


void icy::Element::Reset()
{
    elem_incident_faces.clear();
    elem_incident_faces.reserve(32);
}


void icy::Element::Precompute()
{
    strain_energy_density = 0;
    principal_stress[0] = 0;
    principal_stress[1] = 0;
    principal_stress[2] = 0;
    max_shear_stress = 0;
    hydrostatic_stress = 0;

    computeDm();
    volume_initial = Dm.determinant()/6.;
    if(volume_initial==0) throw std::runtime_error("element's initial area is zero");
    else if(volume_initial < 0)
    {
        std::swap(nds[0],nds[1]);
        computeDm();
        volume_initial = Dm.determinant()/6.;
    }
}



void icy::Element::computeDm()
{
    Dm << (nds[0]->x0 - nds[3]->x0), (nds[1]->x0 - nds[3]->x0), (nds[2]->x0 - nds[3]->x0);
    DmInv = Dm.inverse();
}


int icy::Element::FacetToAbaqusFacetIndex(std::tuple<int,int,int> face)
{

    for(int i=0;i<4;i++)
    {
        int facet[3] = {nds[fi[i][0]]->globId, nds[fi[i][1]]->globId, nds[fi[i][2]]->globId};
        std::sort(std::begin(facet),std::end(facet));
        std::tuple<int,int,int> result(facet[0],facet[1],facet[2]);
        if(face == result) return i;
    }
    spdlog::critical("facet not found in the element");
    for(int i=0;i<4;i++)
    {
        int facet[3] = {nds[fi[i][0]]->globId, nds[fi[i][1]]->globId, nds[fi[i][2]]->globId};
        std::sort(std::begin(facet),std::end(facet));
        std::tuple<int,int,int> result(facet[0],facet[1],facet[2]);
        spdlog::info("{} : {}-{}-{}  vs  {}-{}-{}", i,
                     std::get<0>(face),std::get<1>(face),std::get<2>(face),
                     facet[0],facet[1],facet[2]);
    }

    throw std::runtime_error("facet not found in the element");
}


