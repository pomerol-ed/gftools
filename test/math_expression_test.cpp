#include <numeric>

#include "math_expression.hpp"
#include "container.hpp"
#include "real_grid.hpp"
#include "grid_object.hpp"

#include "gtest/gtest.h"

using namespace gftools;

TEST(math_expr, t1)
{
}

int main(int argc, char **argv) {
    std::cout << "Hi!" << std::endl;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
