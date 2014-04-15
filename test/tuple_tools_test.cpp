#include "defaults.hpp"
#include "tools.hpp"

#include "gtest/gtest.h"

using namespace gftools;
using namespace tuple_tools;

struct tuple_tools_test : public ::testing::Test {
protected:
    virtual void SetUp() {
        tuple1 = std::make_tuple(1,2.,3.0+4.45*I);
    }

    virtual void TearDown() {
    }
public:
    typedef std::tuple<int,double,complex_type> tuple_type;
    tuple_type tuple1;
};



TEST_F(tuple_tools_test, output) {
        std::cout << "tuple:         " << print_tuple(tuple1) << std::endl;
        std::cout << "serialization: " << serialize_tuple(tuple1) << std::endl << std::endl; 
}

TEST_F(tuple_tools_test, compare1) {
    tuple_type tuple2 = tuple1;
    bool check = tools::is_float_equal(tuple1, tuple2);
    EXPECT_EQ(check, 1);
}

TEST_F(tuple_tools_test, streaming) {
    std::stringstream stream1;
    stream1 << serialize_tuple(tuple1);
    tuple_type tuple3 = read_tuple<tuple_type>(stream1);
    bool check = tools::is_float_equal(tuple1, tuple3, 1e-6);
    std::cout <<  std::boolalpha << print_tuple(tuple3) << "==" << print_tuple(tuple1) << " = " << check << std::endl;
    EXPECT_EQ(check, 1);
}

TEST_F(tuple_tools_test, unfold) {
    std::function<double(int,double,complex_type)> f1 = [](int x, double y, complex_type z){return x+y+real(z);};
    std::cout << "function (" << print_tuple(tuple1) << ") = " << unfold_tuple(f1,tuple1) << std::endl;
    struct struct1 { int operator()(int x,double y,complex_type z){return x+y+real(z);}; };
    std::cout << "structure(" << print_tuple(tuple1) << ") = " << unfold_tuple(struct1(),tuple1) << std::endl;
}

TEST(tuple_tools, repeater)
{
    auto tuple1 = repeater<double, 4>::get_tuple(3.0);
    std::cout << "tuple:         " << print_tuple(tuple1) << std::endl;
    auto array1 = repeater<double, 4>::get_array(3.0);
    std::cout << "array:         " << print_array(array1) << std::endl;
}


int main(int argc, char **argv) 
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
