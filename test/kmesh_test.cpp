
#include "kmesh.hpp"

using namespace gftools;
int main()
{
    kmesh k1(32);
    std::cout << k1 << std::endl;

    std::function<real_type(real_type)> F1=[](real_type x) {return cos(x);};
    real_type outF1 = k1.integrate(F1);
    std::cout << "sum_k cos(k) = " << outF1 << std::endl;
    if (!tools::is_float_equal(outF1,0.,1e-8)) return EXIT_FAILURE;

    typedef typename kmesh::point point;
    point p1 = k1.find_nearest(PI);
    std::cout << "PI == " << p1 << std::endl;
    if (!tools::is_float_equal(PI,double(p1),1e-8)) return EXIT_FAILURE;

    kmesh_patch patch1(k1, {0, 2, 4, 8, 14, 16, 30});
    std::cout << "kmesh patch : " << patch1 << std::endl;
    std::vector<int> data_y = {0,1,2};
    std::cout << "y[" << k1[2] << "]=" << patch1.evaluate(data_y,k1[2]) << std::endl;
    if (patch1.evaluate(data_y,k1[2]) != 1) return EXIT_FAILURE;


    std::cout << k1.shift(p1,2.0*PI*31./32.) << std::endl;
    return EXIT_SUCCESS;
}
