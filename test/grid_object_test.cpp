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
    
    auto t1 = gf.get_points({{3, 1, 7}});
    DEBUG(print_tuple(t1));
    auto t2 = gf.get_args({{3, 1, 7}});
    DEBUG(print_tuple(t2));
    auto t3 = gf.get_args(t1);
    EXPECT_EQ(is_float_equal(t3,t2),1);
    EXPECT_EQ(gf.size(), 4*10*8);

    gf_t gf2(gf);
    gf2[0][0][0] = 0.5;
    EXPECT_EQ(gf2.diff(gf), 0.5);
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
}




int main(int argc, char **argv) 
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
    /*
{
    D1.getData()[3]=3.4;
    DEBUG(D1);
    DEBUG(D1(FMatsubara(3,20)));
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> D2(std::make_tuple(n1,n2));
    
    D2.getData()[0][1]=4.0;
    D2.getData()[1][2]=3.1;
    DEBUG(D2);
    D1.savetxt("D1.dat");
    D2.savetxt("D2.dat");

    auto& C1 = D2.getData();
    DEBUG(C1[0]);
    //decltype (C1[0]) x(std::make_tuple(3));
    container<complex_type,1> C2(C1[0]);
    DEBUG(C2[1]);
    DEBUG(n2.evaluate(C2, FMatsubara(1,20)));
    DEBUG(n1.evaluate(C1, FMatsubara(1,10)));
    auto C22 = n1.evaluate(C1, FMatsubara(0,10));
    DEBUG(C22+C2-C22*2.0);
    //DEBUG(n2.getValue(std::forward<Container<complex_type,1>>(n1.getValue(C1, FMatsubara(1,10))), FMatsubara(1,20)));

    DEBUG(D2(FMatsubara(1,10),FMatsubara(2,20)));

    std::function<complex_type(complex_type, complex_type)> f1 = [](const complex_type &a1, const complex_type &a2){return 1.0/(a1+a2);};
    DEBUG(f1(1.0,2.0));
    
    
    D2.fill(f1);
    DEBUG(D2);
    
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid, fmatsubara_grid> D3(std::make_tuple(n1,n2,n2));
    std::function<complex_type(complex_type, complex_type, complex_type)>  f2 = [](const complex_type &a1, const complex_type &a2, const complex_type &a3){return a3/(a1+a2);};
    std::function<complex_type(complex_type, complex_type, complex_type)>  f3 = [](const complex_type &a1, const complex_type &a2, const complex_type &a3){return a3*a2/(a1);};
    D3.fill(f2);
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid, fmatsubara_grid> D4(std::make_tuple(n1,n2,n2));
    D4.fill(f3);
    DEBUG(D4);
    DEBUG(D3);
    DEBUG(D3+D4);

    // Tails test.
    int n_freq = 6;
    real_type beta=10.0;
    fmatsubara_grid fgrid(-n_freq,n_freq,beta);
    real_type t=3.0;
    std::function<complex_type(complex_type)> f11, f21;
    f11 = [t](complex_type w) -> complex_type {return t*t/w;};
    f21 = [t](complex_type w) -> complex_type {return t;};
    grid_object<complex_type,fmatsubara_grid> D21(fgrid);
    grid_object<complex_type,fmatsubara_grid> D31(fgrid);
    D21 = t;
    DEBUG(D21);
    DEBUG(D21(FMatsubara(n_freq+1,beta)));

    D21.fill(f11);
    DEBUG(D21);
    DEBUG(D21(FMatsubara(n_freq+1,beta)));

    
    D21 = 2.0;
    DEBUG(D21(FMatsubara(n_freq+1,beta)));
    D21+=-1.0;
    DEBUG(D21(FMatsubara(n_freq+1,beta)));
    D31 = 3.0;
    DEBUG(D31(FMatsubara(n_freq+1,beta)));

    auto D41 = D21+D31;
    DEBUG(D41);
    DEBUG(D41(FMatsubara(n_freq+1,beta)));

    grid_object<real_type, enum_grid> g1(enum_grid(0,10));
    grid_object<real_type, enum_grid> g2(enum_grid(0,10));
    
    g1.fill(grid_object<real_type, enum_grid>::FunctionType([](int x){return x/2.0;}));
    DEBUG(g1);
    DEBUG(g1(3));
    DEBUG(g1(15));

    g1.savetxt("g1.dat");
    g2.loadtxt("g1.dat");
    return EXIT_SUCCESS;
}
*/
