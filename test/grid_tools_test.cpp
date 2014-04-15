#include "gtest/gtest.h"

#include "matsubara_grid.hpp"
#include "kmesh.hpp"
#include "enum_grid.hpp"

#include "grid_tools.hpp"

using namespace gftools;
using namespace tools;
using namespace tuple_tools;

struct grid_tuple_test : public ::testing::Test {
protected:
    virtual void SetUp() {
        grids_ptr.reset(new grid_tuple(std::make_tuple<enum_grid,fmatsubara_grid,kmesh> (
            enum_grid(0,4),
            fmatsubara_grid(0,12,beta),
            kmesh(8)))
        );
    }

    virtual void TearDown() {
    }
public:
    const double beta = 10.;
    typedef std::tuple<enum_grid,fmatsubara_grid,kmesh> grid_tuple;
    std::unique_ptr<grid_tuple> grids_ptr;
    const grid_tuple& grids(){return *grids_ptr;}
    typedef grid_tuple_traits<grid_tuple> trs;
};

// check construction
TEST_F(grid_tuple_test, Construction) {
    EXPECT_EQ(std::get<0>(*grids_ptr).size(), 4);
    EXPECT_EQ(std::get<1>(*grids_ptr).size(), 12);
    EXPECT_EQ(std::get<2>(*grids_ptr).size(), 8);
}

TEST_F(grid_tuple_test, IndexOp) {
    const auto& g = grids();
    typename trs::point_tuple p = std::make_tuple(std::get<0>(g)[1], std::get<1>(g)[2], std::get<2>(g)[4]);
    std::cout << "p : " << print_tuple(p) << std::endl;
    typename trs::indices i1 = trs::get_indices(p, g);
    std::cout << "indices : " << print_array(i1) << std::endl;
    EXPECT_EQ(i1[0], 1);
    EXPECT_EQ(i1[1], 2);
    EXPECT_EQ(i1[2], 4);
    
    typename trs::indices i2 = {{2,5,3}};
    typename trs::point_tuple p2 = trs::get_points(i2,g);
    typename trs::indices i22 = trs::get_indices(p2, g);
    std::cout << print_array(i2) << "==" << print_array(i22) << std::endl;
    EXPECT_EQ(i2,i22);
    typename trs::arg_tuple a2 = trs::get_args(i2,g);
    typename trs::arg_tuple a22 = p2;
    std::cout << print_tuple(p2) << ":" << print_tuple(a22) << " == " << print_tuple(a2) << std::endl;
    EXPECT_DOUBLE_EQ(std::get<0>(a22),std::get<0>(a2));
    EXPECT_DOUBLE_EQ(std::abs(std::get<1>(a22)),std::abs(std::get<1>(a2)));
    EXPECT_DOUBLE_EQ(std::get<2>(a22),std::get<2>(a2));
    EXPECT_EQ(int(std::get<0>(p2)),std::get<0>(a2));
    EXPECT_DOUBLE_EQ(std::abs(std::complex<double>(std::get<1>(p2))),std::abs(std::get<1>(a2)));
    EXPECT_DOUBLE_EQ(std::get<2>(p2),std::get<2>(a2));

    EXPECT_ANY_THROW( { trs::get_points({{5,11,7}},g); } );
}

TEST_F(grid_tuple_test, Globals) {
    const auto& g = grids();
    auto dims = trs::get_dimensions(g);
    std::cout << print_array(dims) << std::endl;
    EXPECT_EQ(dims[0], 4); 
    EXPECT_EQ(dims[1], 12); 
    EXPECT_EQ(dims[2], 8); 

    size_t s = trs::get_total_size(g);
    std::cout << "size : " << s << std::endl;
    EXPECT_EQ(s,4*12*8);

    }

int main(int argc, char **argv) 
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
