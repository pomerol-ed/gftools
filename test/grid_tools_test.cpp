#include <memory>

#include "gtest/gtest.h"

#include "matsubara_grid.hpp"
#include "kmesh.hpp"
#include "enum_grid.hpp"
#include "real_grid.hpp"

#include "grid_tools.hpp"
// #include "eval_expression.hpp"  // TODO

using namespace gftools;
using namespace tools;
using namespace tuple_tools;

struct grid_tuple_test : public ::testing::Test {
protected:
    virtual void SetUp() {
        grids_ptr.reset(new grid_tuple(std::make_tuple (
            enum_grid(0,4),
            fmatsubara_grid(0,12,beta),
            kmesh(8),
            real_grid(0,6,100,false)))
        );

        grids2_ptr.reset(new grid_tuple2(std::make_tuple (
            enum_grid(0,4),
            kmesh(8)))
        );
    }

    virtual void TearDown() {
    }
public:
    const double beta = 10.;

    typedef std::tuple<enum_grid,fmatsubara_grid,kmesh, real_grid> grid_tuple;
    std::unique_ptr<grid_tuple> grids_ptr;
    const grid_tuple& grids(){return *grids_ptr;}
    typedef grid_tuple_traits<grid_tuple> trs;


    typedef std::tuple<enum_grid,kmesh> grid_tuple2;
    std::unique_ptr<grid_tuple2> grids2_ptr;
    const grid_tuple2& grids2(){return *grids2_ptr;}
    typedef grid_tuple_traits<grid_tuple2> trs2;
};

// check construction
TEST_F(grid_tuple_test, Construction) {
    EXPECT_EQ(std::get<0>(*grids_ptr).size(), 4);
    EXPECT_EQ(std::get<1>(*grids_ptr).size(), 12);
    EXPECT_EQ(std::get<2>(*grids_ptr).size(), 8);
}

TEST_F(grid_tuple_test, IndexOp) {
    const auto& g = grids();
    typename trs::point_tuple p = std::make_tuple(std::get<0>(g)[1], std::get<1>(g)[2], std::get<2>(g)[4], std::get<3>(g)[1]);
    std::cout << "p : " << print_tuple(p) << std::endl;
    typename trs::indices i1 = trs::get_indices(p, g);
    std::cout << "indices : " << print_array(i1) << std::endl;
    EXPECT_EQ(i1[0], 1);
    EXPECT_EQ(i1[1], 2);
    EXPECT_EQ(i1[2], 4);
    
    typename trs::indices i2 = {{2,5,3,1}};
    typename trs::point_tuple p2 = trs::points(i2,g);
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

    EXPECT_ANY_THROW( { trs::points({{5,11,7,1}},g); } );

    typename trs::point_tuple p23 = trs::find_nearest(a2,g);
    std::cout << "p23: " << print_tuple(p23) << std::endl;
    EXPECT_EQ(p23, p2);
}

TEST_F(grid_tuple_test, Globals) {
    const auto& g = grids();
    auto dims = trs::get_dimensions(g);
    std::cout << print_array(dims) << std::endl;
    EXPECT_EQ(dims[0], 4); 
    EXPECT_EQ(dims[1], 12); 
    EXPECT_EQ(dims[2], 8); 
    EXPECT_EQ(dims[3], 100); 

    size_t s = trs::get_total_size(g);
    std::cout << "size : " << s << std::endl;
    EXPECT_EQ(s,4*12*8*100);
}

TEST_F(grid_tuple_test, PointShift) {
    const auto& g2 = grids2(); 
    trs2::point_tuple p1 = trs2::find_nearest(std::make_tuple(2,PI),g2);
    std::cout << "p1" << print_tuple(p1) << std::endl;

    trs2::arg_tuple shift1 = std::make_tuple(1,PI);
    trs2::point_tuple shift2 = trs2::find_nearest(shift1,g2);

    trs2::point_tuple r1 = trs2::find_nearest(std::make_tuple(3,0.),g2);
    
    trs2::point_tuple p1s1 = trs2::shift(p1,shift1,g2);
    std::cout << print_tuple(p1) << "+" << print_tuple(shift1) << " = " << print_tuple(p1s1) << std::endl;
    EXPECT_EQ(p1s1, r1);

    trs2::point_tuple p1s2 = trs2::shift(p1,shift2,g2);
    std::cout << print_tuple(p1) << "+" << print_tuple(shift2) << " = " << print_tuple(p1s2) << std::endl;
    EXPECT_EQ(p1s2, r1);
}

TEST_F(grid_tuple_test, ArgShift) {
    const auto& g = grids(); 
    trs::point_tuple p1 = trs::find_nearest(std::make_tuple(2,FMatsubara(2,beta),PI, 3. ),g);

    trs::arg_tuple shift1 = std::make_tuple(1,BMatsubara(0,beta),5*PI/4.,4.);
    ASSERT_ANY_THROW(trs::shift(p1,shift1,g)); // last index out of bounds

    shift1 = std::make_tuple(1,BMatsubara(3,beta),5*PI/4.,0.48);
    trs::arg_tuple a1 = std::make_tuple(3,FMatsubara(5,beta),PI/4.,3.48);
    trs::point_tuple r1 = trs::find_nearest(a1,g);

    trs::point_tuple p1s1 = trs::shift(p1,shift1,g);
    std::cout << print_tuple(p1) << "+" << print_tuple(shift1) << " = " << print_tuple(p1s1) << " == " << print_tuple(r1) << std::endl;
    EXPECT_EQ(p1s1, r1);

    trs::arg_tuple a2 = trs::shift(p1,shift1,g);
    DEBUG(print_tuple(a1) << " == " << print_tuple(a2));
    EXPECT_EQ(is_float_equal(a1,a2),1);

    trs::point_tuple p1s2 = trs::find_nearest(a2,g);
    EXPECT_EQ(p1s2, r1);
}

TEST_F(grid_tuple_test, ArgShiftFail) {
    const auto& g = grids(); 
    trs::point_tuple p1 = trs::find_nearest(std::make_tuple(2,FMatsubara(2,beta),PI, 3. ),g);

    double rstep = (std::get<3>(g)[1] - std::get<3>(g)[0]);

    // last index out of bounds
    ASSERT_ANY_THROW(trs::shift(p1,std::make_tuple(1,BMatsubara(0,beta),5*PI/4.,rstep*100),g)); 
    // can't shift fermionic Matsubara over another fermionic Matsubara and get a Fermionic Matsubara
    ASSERT_ANY_THROW(trs::shift(p1,std::make_tuple(1,FMatsubara(0,beta),5*PI/4.,rstep),g)); 
    // first index out of bounds
    ASSERT_ANY_THROW(trs::shift(p1,std::make_tuple(4,BMatsubara(0,beta),5*PI/4.,rstep),g)); 

    // grid point mismatch - one can only shift to the exact point on the grid. Otherwise use find_nearest(arg_tuple)
    ASSERT_ANY_THROW(trs::shift(p1,std::make_tuple(0,BMatsubara(0,beta),5*PI/8.,rstep),g)); 
    // grid point mismatch - one can only shift to the exact point on the grid. Otherwise use find_nearest(arg_tuple)
    ASSERT_ANY_THROW(trs::shift(p1,std::make_tuple(0,BMatsubara(0,beta),5*PI/4.,rstep/2),g)); 

    trs::arg_tuple a1 = trs::shift(trs::get_args(p1,g), std::make_tuple(0,0.,0.,rstep / 3.),g);
    EXPECT_EQ(trs::find_nearest(a1,g),p1);
}

/* // TODO
TEST(grids, EvaluateExprTemp1) {
    
    int size = 10;
    std::vector<double> d(size);
    for (int i=0; i<size; i++) d[i]=double(i)/(size);
    kmesh k1(size);

    grid_eval_expr<decltype(d), std::tuple<kmesh>> expr1 (d,std::forward_as_tuple(k1)); 
    //EXPECT_EQ(d[3], expr1[3]);
    DEBUG((std::is_same<std::vector<double>, std::vector<double>>::value));

    std::vector<std::vector<double>> d2(size);
    for (int i=0; i<size; i++) { d2[i].resize(size); for (int j=0; j<size; j++) d2[i][j] = double(i+j)/size/size; }
    
    //grid_eval_expr<decltype(d2), std::tuple<kmesh,kmesh>> expr2 (d2,std::forward_as_tuple(k1,k1)); 
    //DEBUG(d2[2][3]);
    //DEBUG(expr2[2][3]);
    
    }
*/


int main(int argc, char **argv) 
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
