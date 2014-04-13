
#include "enum_grid.hpp"

using namespace gftools;

int main()
{
    enum_grid g1(3,8);
    INFO(g1);
    enum_grid::point p1 = g1[0];
    std::cout << p1 << std::endl;
    return EXIT_SUCCESS;
}
