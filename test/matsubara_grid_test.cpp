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

    complex_type freq = FMatsubara(10,20);
    bool success = tools::is_float_equal(freq, n2.find_nearest(freq).val_);
    std::cout << freq << " == " << n2.find_nearest(freq) << " = " << std::boolalpha << success << std::endl;
    if (!success) return EXIT_FAILURE;

    freq = FMatsubara(20,20);
    success = tools::is_float_equal(freq, n2.find_nearest(freq).val_);
    std::cout << freq << " == " << n2.find_nearest(freq) << " = " << std::boolalpha << success << std::endl;
    if (!success) return EXIT_FAILURE;

    freq = FMatsubara(32,20);
    success = (n2.size()-1 == n2.find_nearest(freq).index_);
    if (!success) return EXIT_FAILURE;

    INFO(n2.find_nearest(FMatsubara(35,20)));

    std::vector<int> x(n2.size());
    for (size_t i=0; i<x.size(); ++i) x[i]=i*i;
    INFO(n2.evaluate(x,FMatsubara(10,20)));
    if (n2.evaluate(x,FMatsubara(10,20))!=100) return EXIT_FAILURE;

    bmatsubara_grid bgrid1(-20,20,10);

    return EXIT_SUCCESS;
}
