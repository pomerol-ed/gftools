#include <numeric>

#include "Defaults.hpp"
#include "MatsubaraGrid.hpp"
#include "KMesh.hpp"

#include <iostream>
#include <ctime>
#include <array>

using namespace GFTools;

int main()
{
    FMatsubaraGrid n1(-100,100,10);
    FMatsubaraGrid n2(0,32,20);

    KMesh k1(32);
    std::function<RealType(RealType)> F1=[](RealType x) {return cos(x);};
    RealType outF1 = k1.integrate (F1);
    INFO(outF1);
    if (std::abs(outF1)>1e-8) return EXIT_FAILURE;

    std::function<ComplexType(ComplexType)> F2=[](ComplexType x) {return x;};
    ComplexType out = n1.integrate(F2);
    if (std::abs(out)>1e-8) return EXIT_FAILURE;
    INFO(out);
    DEBUG(n2);

    auto a1 = std::make_tuple(n1,n2);
    
    std::vector<int> x(n2.getSize());
    for (int i=0; i<x.size(); ++i) x[i]=i*i;

    INFO(std::get<1>(n2.find(FMatsubara(10,20))));
    if (std::get<1>(n2.find(FMatsubara(10,20)))!=10) return EXIT_FAILURE;
    INFO(std::get<1>(n2.find(FMatsubara(20,20))));
    if (std::get<1>(n2.find(FMatsubara(20,20)))!=20) return EXIT_FAILURE;
    INFO(std::get<1>(n2.find(FMatsubara(32,20))));
    if (std::get<1>(n2.find(FMatsubara(32,20)))!=32) return EXIT_FAILURE;
    INFO(std::get<1>(n2.find(FMatsubara(35,20))));

    INFO(n2.getValue(x,FMatsubara(10,20)));
    if (n2.getValue(x,FMatsubara(10,20))!=100) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
