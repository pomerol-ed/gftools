#include <numeric>

#include "matsubara_grid.hpp"

using namespace gftools;

int main()
{

    fmatsubara_grid n1(-100,100,10);
    fmatsubara_grid n2(0,32,20);

    typename fmatsubara_grid::point ppp = n1[0];
    std::cout << ppp << std::endl;    

    std::function<complex_type(complex_type)> F2=[](complex_type x) {return x;};
    complex_type out = n1.integrate(F2);
    if (std::abs(out)>1e-8) return EXIT_FAILURE;
    INFO(out);
    DEBUG(n2);

    auto a1 = std::make_tuple(n1,n2);
    
    std::vector<int> x(n2.size());
    for (size_t i=0; i<x.size(); ++i) x[i]=i*i;

    INFO(std::get<1>(n2.find(FMatsubara(10,20))));
    if (std::get<1>(n2.find(FMatsubara(10,20)))!=10) return EXIT_FAILURE;
    INFO(std::get<1>(n2.find(FMatsubara(20,20))));
    if (std::get<1>(n2.find(FMatsubara(20,20)))!=20) return EXIT_FAILURE;
    INFO(std::get<1>(n2.find(FMatsubara(32,20))));
    if (std::get<1>(n2.find(FMatsubara(32,20)))!=32) return EXIT_FAILURE;
    INFO(std::get<1>(n2.find(FMatsubara(35,20))));

    INFO(n2.evaluate(x,FMatsubara(10,20)));
    if (n2.evaluate(x,FMatsubara(10,20))!=100) return EXIT_FAILURE;



    return EXIT_SUCCESS;
}
