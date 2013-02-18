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
    Log.setDebugging(true);
    size_t n_freq = 5;
    RealType beta = 10;
    FMatsubaraGrid gridF(0, n_freq, beta);
    BMatsubaraGrid gridB(0, n_freq, beta);

    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1;
    RealType t = 0.5;
    f1 = [t](ComplexType w) -> ComplexType {return t*t/w;};
    Delta.fill(f1);
    
    auto D2 = Delta.shift(gridB[1]);
    if (D2[0]!=Delta[1]) return EXIT_FAILURE;

    KMesh kGrid(8);
    DEBUG(kGrid.shift(PI*5.0/4,PI) - PI/4.);
    DEBUG(kGrid[7] << "|" << 7.0/8*2*PI);
    DEBUG(kGrid.shift(7./4*PI,PI) << "|" << 3/4.*PI);
    DEBUG(kGrid.shift(kGrid[7],PI) << "|" << 3/4.*PI << "|" << kGrid[3]);
    bool success = (kGrid.shift(kGrid[7],PI)==kGrid[3]);
    if (!success) { ERROR("Point mismatch after shift."); return EXIT_FAILURE; };

    GridObject<RealType,KMesh,KMesh> e_k(std::forward_as_tuple(kGrid,kGrid));
    auto ekf = [](RealType kx, RealType ky){return 2.0*(cos(kx)+cos(ky));};
    e_k.fill(ekf);
    DEBUG(e_k);
    DEBUG(e_k.shift(PI,PI));

    bool b1 = is_equal((e_k.shift(PI,PI)+e_k).sum(),0);
    if (!b1) return EXIT_FAILURE;

    kGrid = KMesh(32);
    DEBUG(kGrid.shift(kGrid[31],PI) << "==" << 15./16*PI);
    DEBUG(kGrid.shift(kGrid[31],-5*PI) << "==" << 15./16*PI);
    
    return EXIT_SUCCESS;
}
