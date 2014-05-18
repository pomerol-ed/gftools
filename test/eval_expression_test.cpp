#include <numeric>

#include "eval_expression.hpp"
#include "container.hpp"
#include "real_grid.hpp"
#include "grid_object.hpp"

#include "gtest/gtest.h"

using namespace gftools;

TEST(expression1d, vector) {
    std::vector<double> y(10);
    std::iota(y.begin(), y.end(), 0.3);
    eval_expression<std::vector<double>, double, 1> e(y);
    std::cout <<  y[2] << "==" << e[2] << std::endl;
    for (int i=0; i<y.size(); ++i) { 
        ASSERT_EQ(y[i], e[i]);
        };
}

TEST(expression1d, container) {
    container<double, 1> y(10);
    std::iota(&y[0], &y[9], 0.3);
    eval_expression<container<double, 1>, double, 1> e(y);
    std::cout <<  y[2] << "==" << e[2] << std::endl;
    for (int i=0; i<y.size(); ++i) { 
        ASSERT_EQ(y[i], e[i]);
        };
}

TEST(expression2d, vector) {
    std::array<size_t, 2> dims = {{10, 4}};
    typedef std::vector<std::vector<double>> data2d;
    data2d y(dims[0]);
    for (auto& x : y) { x.resize(dims[1]); std::iota(x.begin(), x.end(), 1.2); };
    eval_expression<data2d, double, 2> e(y);
    std::cout << "original   : " << y[2][3] << std::endl;
    std::cout << "expression : " << e[2][3] << std::endl;
    DEBUG(y[2][3] + y[1][2]);
    DEBUG(e[2][3] + e[1][2]);
    for (int i1 = 0; i1 < dims[0]; i1++) 
        for (int i2 = 0; i2 < dims[1]; i2++) 
            ASSERT_EQ(y[i1][i2], e[i1][i2]);
}

TEST(expression2d, container) {
    typedef container<double, 2> data2d;
    std::array<size_t, 2> dims = {{10, 4}};
    data2d y(dims);
    std::iota(&y[0][0], &y[dims[0]-1][dims[1]-1], 0.2);
    eval_expression<data2d, double, 2> e(y);
    std::cout << "original   : " << y[2][3] << std::endl;
    std::cout << "expression : " << e[2][3] << std::endl;
    for (int i1 = 0; i1 < dims[0]; i1++) 
        for (int i2 = 0; i2 < dims[1]; i2++) 
            ASSERT_EQ(y[i1][i2], e[i1][i2]);
    };

TEST(expression2d, grid_object) {
    typedef grid_object<double, real_grid, real_grid> data2d;
    std::array<size_t, 2> dims = {{10, 4}};
    data2d y(real_grid(0,1.,dims[0]), real_grid(0,1.,dims[1]));
    y.fill(typename data2d::function_type([](double x, double y){return sin(x+y);}));
    eval_expression<data2d, double, 2> e(y);
    std::cout << "original   : " << y[2][3] << std::endl;
    std::cout << "expression : " << e[2][3] << std::endl;
    for (int i1 = 0; i1 < dims[0]; i1++) 
        for (int i2 = 0; i2 < dims[1]; i2++) 
            ASSERT_EQ(y[i1][i2], e[i1][i2]);
    };


TEST(expression3d, vector) {
    typedef std::vector<std::vector<std::vector<double>>> data3d;
    std::array<size_t, 3> dims = {{10, 4, 6}};
    data3d y(dims[0]);
    for (auto& x : y) { x.resize(dims[1]); for (auto& z: x) { z.resize(dims[2]); std::iota(z.begin(), z.end(), 1.2); };};
    eval_expression<data3d, double, 3> e(y);
    std::cout << "original   : " <<  y[2][3][5] << std::endl;
    std::cout << "expression : " <<  e[2][3][5] << std::endl;
    for (int i1 = 0; i1 < dims[0]; i1++) 
        for (int i2 = 0; i2 < dims[1]; i2++) 
            for (int i3 = 0; i3 < dims[2]; i3++) 
                ASSERT_EQ(y[i1][i2][i3], e[i1][i2][i3]);
}

TEST(expression3d, container) {
    typedef container<double, 3> data3d;
    std::array<size_t, 3> dims = {{10, 4, 6}};
    data3d y(dims);
    std::iota(&y[0][0][0], &y[dims[0]-1][dims[1]-1][dims[2]-1], 0.1);
    eval_expression<data3d, double, 3> e(y);
    std::cout << "original   : " <<  y[2][3][5] << std::endl;
    std::cout << "expression : " <<  e[2][3][5] << std::endl;
    for (int i1 = 0; i1 < dims[0]; i1++) 
        for (int i2 = 0; i2 < dims[1]; i2++) 
            for (int i3 = 0; i3 < dims[2]; i3++) 
                ASSERT_EQ(y[i1][i2][i3], e[i1][i2][i3]);
}



int main(int argc, char **argv) {
    std::cout << "Hi!" << std::endl;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
