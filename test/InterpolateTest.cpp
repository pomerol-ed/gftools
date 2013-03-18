#include <numeric>

#include "Defaults.hpp"
#include "MatsubaraGrid.hpp"
#include "KMesh.hpp"
#include "GridObject.hpp"

#include <iostream>
#include <ctime>
#include <array>

using namespace GFTools;
typedef GridObject<ComplexType,FMatsubaraGrid> GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}


int main()
{
    RealType beta=10;
    size_t n_freq = 5;
    size_t n_freq2 = 10;
    FMatsubaraGrid gridF(0, n_freq, beta);
    FMatsubaraGrid gridF2(0, n_freq2, beta);

    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1;
    RealType t = 0.5;
    f1 = [t](ComplexType w) -> ComplexType {return t*t/w;};
    Delta.fill(f1);
    GF Delta2(gridF2);
    GF Delta3(gridF2);
    Delta3.fill(f1);
    
    Delta2.copyInterpolate(Delta);
    DEBUG(Delta);
    DEBUG(Delta2);
    DEBUG(Delta3);
    if (!is_equal(Delta2.diff(Delta3),0.0)) return EXIT_FAILURE;

    KMesh kGrid(2);
    GridObject<ComplexType, FMatsubaraGrid, KMesh, KMesh> G1(std::make_tuple(gridF,kGrid,kGrid));
    GridObject<ComplexType, FMatsubaraGrid, KMesh, KMesh> G2(std::make_tuple(gridF2,kGrid,kGrid));
    GridObject<ComplexType, FMatsubaraGrid, KMesh, KMesh> G3(std::make_tuple(gridF2,kGrid,kGrid));

    std::function<ComplexType(ComplexType,RealType,RealType)> f2 = [t](ComplexType w, RealType kx, RealType ky){return 1.0/(w+2*t*cos(kx)+2*t*cos(ky));};
    G1.fill(f2);
    G3.fill(f2);
    G2.copyInterpolate(G1);

    DEBUG(G1); 
    DEBUG(G2); 
    DEBUG(G3); 

    DEBUG(G2.diff(G3));
    if (!is_equal(G2.diff(G3),0.0)) return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}
