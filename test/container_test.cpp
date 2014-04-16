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

    EXPECT_DOUBLE_EQ((*cd2).diff(cd2->conj().conj()), 0);
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
    std::cout << "Hi!" << std::endl;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
