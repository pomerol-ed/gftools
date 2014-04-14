#include <numeric>
#include <iostream>
#include <ctime>
#include <array>
#include <functional>

#include "defaults.hpp"
#include "container.hpp"

#include "gtest/gtest.h"

using namespace gftools;

TEST(ContainerTest, ConstructionFromDims) {
    container<double,1> d1(5);
    container<double,2> d2(5,4);
    container<double,3> d3(5,4,2);

    container<std::complex<double>,1> cd1(5);
    container<std::complex<double>,2> cd2(5,4);
    container<std::complex<double>,3> cd3(5,4,2);

    container<int,1> i1(5);
    container<int,2> i2(5,4);
    container<int,3> i3(5,4,2);
}

struct cont_test : public ::testing::Test {
protected:
    virtual void SetUp() {
        d1.reset(new container<double,1>(shape1));
        d2.reset(new container<double,2>(shape2));
        d3.reset(new container<double,3>(shape3));

        cd1.reset(new container<std::complex<double>,1>(shape1));
        cd2.reset(new container<std::complex<double>,2>(shape2));
        cd3.reset(new container<std::complex<double>,3>(shape3));
    }

    virtual void TearDown() {
    }
public:
    std::array<size_t,1> shape1 {{5}};
    std::array<size_t,2> shape2 {{3,3}};
    std::array<size_t,3> shape3 {{2,3,3}};

    std::unique_ptr<container<double,1>> d1;
    std::unique_ptr<container<double,2>> d2;
    std::unique_ptr<container<double,3>> d3;

    std::unique_ptr<container<std::complex<double>,1>> cd1;
    std::unique_ptr<container<std::complex<double>,2>> cd2;
    std::unique_ptr<container<std::complex<double>,3>> cd3;
};

// check construction
TEST_F(cont_test, ConstructionFromShape) {
    EXPECT_EQ(d1->shape(), shape1);
    EXPECT_EQ(d2->shape(), shape2);
    EXPECT_EQ(d3->shape(), shape3);
    EXPECT_EQ(cd2->shape(), shape2);
}


TEST_F(cont_test, Assignment) {
    (*d2) = 2.0;
    EXPECT_DOUBLE_EQ((*d2).sum(), std::accumulate(shape2.begin(), shape2.end(), 1, std::multiplies<int>())*2 );

    (*d3) = 1.0;
    auto a1 = (*d3)[0] - (*d2);

    (*d2)+=a1;

    container<double,2> t1(shape2);
    t1 = -1.;
    EXPECT_EQ(a1.sum(), t1.sum());
}

TEST_F(cont_test, Math) {
    (*d2) = 2.0;
    (*d3) = 1.0;
    auto a1 = (*d3)[0] - (*d2);
    auto a2 = 3.0*(a1+0.5)*2+(*d2);

    EXPECT_DOUBLE_EQ(a2[0][0], -1.);

    (*cd2) = 3.0+ 4.0*I;
    (*cd2)*=cd2->conj();
    EXPECT_DOUBLE_EQ(std::abs((*cd2)[0][0]), 25); 
}

TEST_F(cont_test, Iterators) {
    (*d2) = 0.0;
    (*d2)[0][1] = 3.0;
    (*d2)[0][2] = 4.0;

    auto it1 = d2->begin()->begin();
    std::advance(it1,1);
    EXPECT_DOUBLE_EQ(*it1, 3.0);

    (*d1) = 1.0;
}

TEST_F(cont_test, Flatten) {
    auto f = (*d2).flatten();
    EXPECT_EQ(f.size(), (*d2).size());
    EXPECT_EQ(f.shape()[0], (*d2).size());
}


int main(int argc, char **argv) {

    DEBUG(bool(std::is_convertible<std::tuple<int,int>, typename tuple_tools::repeater<int,2>::tuple_type >::value));

    std::cout << "Hi!" << std::endl;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
    
/*
    std::array<size_t,1> Ar1 {{5}};
    std::array<size_t,2> Ar2 {{3,3}};
    std::array<size_t,3> Ar3 {{2,3,3}};


    INFO("Container")

    container<complex_type,2> A1(Ar2), A2(Ar2), A3(Ar2);
    A1 = 1.;
    A2 = 2.;
    DEBUG(A1);
    A1+=A2;
    DEBUG(A1);
    A3 = 3.0*(A1+0.5)*2+A2;
    DEBUG(A3);
    DEBUG(A1);
    auto a01 = 3*(A1[1]+A1[0])*2;
    DEBUG(a01);
    DEBUG(A1);

    container<complex_type,1> E1(A1[1]), E2(A1[2]);
    E1*=3;
    DEBUG(E1);
    DEBUG(A1);
    

    container<complex_type,2> B(Ar2);
    container<complex_type,2> C(Ar2);
    container<complex_type,3> D(Ar3);
    
    B[0][2]=3.0;
    C[0][2]=-2.0;

    std::array<size_t, 2> coord1 {{0,1}};
    DEBUG(B.data_(coord1));

    
    container<complex_type,2>::iterator t1;
    auto it2 = B.begin();
    auto it3 = B[0].begin();
    
    DEBUG(B[0][2]);
    DEBUG(B[0][0]);

    //DEBUG(B.data_.shape()[0]);
    //DEBUG(B[0].data_.shape()[0]);

    DEBUG(B);
    DEBUG(B[0]);

    D[0][0][0]=-1.0;
    D[0][0][1]=1.0;
    DEBUG(D);
    DEBUG(D.data_.num_elements());

    DEBUG(C);
    C=B;
    DEBUG(C);

    C[0][0]=1.0;
    C[0][2]=-5;
    DEBUG(C);
    C+=B;
    DEBUG(C);
    C+=1.0;
    DEBUG(C);
    DEBUG(C+B);
    container<complex_type,3> E(D);

    DEBUG(B);
    B+=2.0;
    DEBUG(B);
    DEBUG(B[0]+C[0]);
    DEBUG((B+C)[0]);
    DEBUG((B-C)[0]);

    B*=4.;
    DEBUG(B);
    DEBUG(B[0]*2.0);
    DEBUG(D+2.0);
    DEBUG(D[0]+B);
    //DEBUG((B[0]*3.0));
    E*=(-1);
    DEBUG(D);
    DEBUG(E);
    DEBUG(D+E);
    DEBUG(D*5+2.0);

    DEBUG(D.conj());


    
    INFO("Matrix test");
    INFO("2-dim");
    std::array<size_t,2> dim1 {{4,5}};
    container<double,2> Vals (dim1);
    Vals[1][0] = 1.0;
    Vals[0][3] = -0.4;
    Vals[2][1] = 1.5;
    INFO(Vals);
    INFO(Vals.as_matrix());
    DEBUG("!");
    container<double,2> Vals_2(Vals.as_matrix());
    INFO(Vals_2);

    decltype(Vals) Vals_3(dim1);
    Vals_3 = Vals.as_matrix();
    INFO(Vals_3);

    if ((Vals-Vals_2).sum()!=0 || (Vals-Vals_3).sum()!=0) return EXIT_FAILURE;
    INFO("1-dim test");
    container<double,1> Vals_1d(std::array<size_t,1>({{5}}));
    Vals_1d[2]=1.6;
    Vals_1d[4]=3;
    INFO(Vals_1d);
    INFO(Vals_1d.as_vector());
    INFO(Vals_1d.as_diagonal_matrix());

    INFO("View test");
    container<double, 3> C1(std::array<size_t,3>({{2,1,3}}));
    C1[0][0][1] = 1.0; C1[0][0][2] = 2.0; C1[1][0][0] = 3.0; C1[1][0][2] = 4.0;
    INFO_NONEWLINE("Data in memory: "); for (size_t i = 0; i<C1.data_.num_elements(); ++i) { INFO_NONEWLINE(*(C1.data_.origin()+i) << " "); }; INFO("");
    INFO("container in : " << C1);
    typedef boost::multi_array<double, 3>::array_view<3>::type myview_3_t;
    typedef boost::multi_array_types::index_range range;
    myview_3_t boost_view1 = C1.data_[boost::indices[range(0,1)][range(0,1)][range(0,2)]]; 
    container_view<double, 3, 3> C1_v2(myview_3_t(C1.data_[boost::indices[range().stride(2)][range()][range().stride(1)]]));
    container_view<double, 3, 3> C1_v3(C1_v2);
    INFO(C1_v2);
    DEBUG(C1_v2*C1_v3);
    DEBUG(C1);
    C1_v2 *= C1_v3;
    DEBUG(C1_v2);
    DEBUG(C1);
    container_view<double, 3> C1_v(C1.data_[boost::indices[range()][range()][range().stride(3)]]);
    INFO("Container view : " << C1_v);

*/
    return EXIT_SUCCESS;
}
