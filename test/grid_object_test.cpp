#include <numeric>

#include "gtest/gtest.h"

#include "defaults.hpp"
#include "matsubara_grid.hpp"
#include "container.hpp"
#include "grid_object.hpp"
#include "kmesh.hpp"
#include "enum_grid.hpp"

using namespace gftools;
using namespace tools;
using namespace tuple_tools;

TEST(GridObject1, Construction1dAndRef) {

    fmatsubara_grid grid1(0,5,10);
    typedef grid_object<complex_type,fmatsubara_grid> gf;

    grid_object_base<container<complex_type,1>,fmatsubara_grid> D11(grid1);
    gf D1(grid1);
    D1.data()[1]=1.1; D1.data()[3]=4.7;
    grid_object_base<container_ref<complex_type,1>,fmatsubara_grid> D12(std::forward_as_tuple(grid1),D1.data());

    EXPECT_DOUBLE_EQ(std::abs(D1.data().sum()), std::abs(D12.data().sum()));

    D1.data()[2] = 2.5;
    EXPECT_DOUBLE_EQ(std::abs(D1.data().sum()), std::abs(D12.data().sum()));

    // from move constructor of container
    gf D3(std::forward_as_tuple(grid1), container<complex_type, 1>(5));
    // provide mismatch in size of grid and container
    ASSERT_ANY_THROW(gf(std::forward_as_tuple(grid1), container<complex_type, 1>(3)));

    // move constructor
    gf(std::forward<gf>(gf(std::forward_as_tuple(grid1), container<complex_type, 1>(5))));
}

TEST(GridObject1, Construction1dFromContainer) {
    fmatsubara_grid grid1(0,5,10);
    typedef grid_object<complex_type,fmatsubara_grid> gf;
    container<complex_type, 1> data1 (5);
    data1[1] = 1.98; data1[2] = 3.17;
    gf D2(std::forward_as_tuple(grid1), data1);

    grid_object_base<container_ref<complex_type,1>,fmatsubara_grid> D12(std::forward_as_tuple(grid1),data1);
    grid_object_base<container_ref<complex_type,1>,fmatsubara_grid> D22(std::forward_as_tuple(grid1),D2.data());

    D2.data()[1] = 2.0;
    data1[1] = 4.0;
  
    // now D12 = data1 and D22 = D2, but not D12 = D2 
    EXPECT_EQ(is_float_equal(D2.data().sum(), data1.sum()),0); 
    EXPECT_DOUBLE_EQ(std::abs(D12.data().sum()), std::abs(data1.sum())); 
    EXPECT_EQ(is_float_equal(D12.data().sum(), D2.data().sum()),0); 
    EXPECT_DOUBLE_EQ(std::abs(D22.data().sum()), std::abs(D2.data().sum())); 
    EXPECT_EQ(is_float_equal(D12.data().sum(), D2.data().sum()),0); 
}

TEST(GridObject1, Construction2d) {
    fmatsubara_grid grid1(0,5,10);
    typedef grid_object<double,fmatsubara_grid, fmatsubara_grid> gf2;
    gf2 D1(std::forward_as_tuple(grid1,grid1));
    gf2 D12(grid1,grid1);
    D1.data()[0][3] = 3.0;
    grid_object_ref<double,fmatsubara_grid,fmatsubara_grid> D3 (std::forward_as_tuple(grid1,grid1), D1.data());
    ASSERT_DOUBLE_EQ(D3.data()[0][3], 3.0);

    D12 = D1;
    EXPECT_DOUBLE_EQ(D12.diff(D1),0);

    fmatsubara_grid grid2(0,4,10);
    // can't assign object of different size
    ASSERT_ANY_THROW(D12 = gf2(grid1,grid2));
}

struct gridobject_test1 : public ::testing::Test {
protected:
    virtual void SetUp() {
        gf_ptr.reset(new gf_t(std::forward_as_tuple(enum_grid(0,4),fmatsubara_grid(0,10,beta),kmesh(8))));
    }

    virtual void TearDown() {
    }
public:
    typedef grid_object<double,enum_grid,fmatsubara_grid,kmesh> gf_t;
    typedef typename gf_t::point_tuple point_tuple;
    typedef typename gf_t::arg_tuple arg_tuple;
    double beta = 10.;
    std::unique_ptr<gf_t> gf_ptr;
};

TEST_F(gridobject_test1, Reductions)
{
    gf_t& gf = *gf_ptr;
    gf[0][2][3] = 1.0;
    gf[1][3][2] = 4.0;
    gf[3][1][1] = 5.0;

    gf.data()[0][1][2] = 2.0;

    EXPECT_DOUBLE_EQ(gf.sum(),12.);

    gf.tail() = [](int a, complex_type b, double c){return 1.0;}; 
    EXPECT_DOUBLE_EQ(gf.tail()(1.,I,3.),1.0);
    
    auto t1 = gf.points({{3, 1, 7}});
    DEBUG(print_tuple(t1));
    auto t2 = gf.get_args({{3, 1, 7}});
    DEBUG(print_tuple(t2));
    auto t3 = gf.get_args(t1);
    EXPECT_EQ(is_float_equal(t3,t2),1);
    EXPECT_EQ(gf.size(), 4*10*8);

    gf_t gf2(gf);
    gf2[0][0][0] = 0.5;
    EXPECT_DOUBLE_EQ(gf2.diff(gf, false), 0.5);
}

TEST_F(gridobject_test1, IO)
{
    gf_t& gf = *gf_ptr;
    gf[0][2][3] = 1.0;
    gf[1][3][2] = 4.0;
    gf[3][1][1] = 5.0;

    gf.savetxt("g1.dat");
    gf_t g2(gf.grids());
    g2.loadtxt("g1.dat");
    EXPECT_DOUBLE_EQ(g2.diff(gf),0.0);

    auto g3 = loadtxt<gf_t>("g1.dat");
    std::cout << "(g3 == g) = " << is_float_equal(g3.diff(gf), 0.0) << std::endl;
    g3.savetxt("g2.dat");
    EXPECT_DOUBLE_EQ(g3.diff(gf),0.0);
}


TEST_F(gridobject_test1, fill)
{
    gf_t& gf = *gf_ptr;
    gf_t gf2(gf.grids());

    typename gf_t::function_type f1  = []  (int_wrap_enumerate_grid i, complex_type p, double k)->double{return int(i)*std::abs(p) + k;};
    gf.fill(f1);

    typename gf_t::point_function_type f2 = [f1](typename enum_grid::point i, typename fmatsubara_grid::point p, typename kmesh::point k){return f1(i,p,k);};
    gf2.fill(f2);
    EXPECT_DOUBLE_EQ(gf2.diff(gf),0.0);

    std::function<double(point_tuple)> f3 = [&](point_tuple in){return unfold_tuple(f2,in);};
    gf2.fill(f3);
    EXPECT_DOUBLE_EQ(gf2.diff(gf),0.0);

    gf2 = 0.0;
    //gf2.copy_interpolate(gf);
}

TEST_F(gridobject_test1, shift)
{
    gf_t& gf = *gf_ptr;
    gf = 1.0;

    gf_t gf2 = gf.shift(0,0.*I,0.);
    EXPECT_DOUBLE_EQ(gf2.diff(gf),0);

    gf.fill([](int i, complex_type p, double k){return cos(k/2.);});
    gf2 = gf.shift(0,0,PI);
    
    EXPECT_DOUBLE_EQ((gf2*gf2+gf*gf).sum(),gf2.size());
}




int main(int argc, char **argv) 
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

