
#include "real_grid.hpp"

using namespace gftools;
int main()
{
    real_grid grid1(-5,5,20);
    std::cout << "grid1 : " << grid1 << std::endl;
    std::cout << "grid1 is uniform : " << std::boolalpha << grid1.check_uniform_() << std::endl;
    
    std::function<real_type(int)> f1;
    f1 = [](int n){return std::pow(3,n);};
    real_grid grid2(0,10,f1);
    std::cout << "grid2 : " << grid2 << std::endl;
    std::cout << "grid2 is uniform : " << std::boolalpha << grid2.check_uniform_() << std::endl;
    
    std::vector<real_type> v1;
    for (real_type a=-2;a<=4;a+=0.1) v1.push_back(a);
    real_grid grid3(v1);
    std::cout << "grid3 : " << grid3 << std::endl;
    std::cout << "grid3 is uniform : " << std::boolalpha << grid3.check_uniform_() << std::endl;
    
    std::function<real_type(real_type)> sinF;
    sinF = [](real_type x){return sin(x);};
    auto int_val = (grid1.integrate(sinF));
    std::cout << "int_grid1 sin(x) dx = " << int_val << std::endl;
    if (!tools::is_float_equal(int_val, 0., 1e-8)) return EXIT_FAILURE;

    auto p1 = grid3.find_nearest(2.);
    std::cout << "2 = " << p1 << std::endl;
    
    return EXIT_SUCCESS;
}
