#include "defaults.hpp"
#include "tools.hpp"

using namespace gftools;
using namespace tuple_tools;

int main()
{
    {
        typedef std::tuple<int,double,complex_type> tuple_type;
        tuple_type tuple1 = std::make_tuple(1,2.,3.0+4.45*I);
        std::cout << "tuple:         " << print_tuple(tuple1) << std::endl;
        std::cout << "serialization: " << serialize_tuple(tuple1) << std::endl << std::endl; 

        std::cout << "TEST : tuple compare" << std::endl;
        tuple_type tuple2 = tuple1;
        bool check = tools::is_float_equal(tuple1, tuple2);
        std::cout << std::boolalpha << print_tuple(tuple1) << "==" << print_tuple(tuple2) << " = " << check << std::endl; 
        if (!check) return EXIT_FAILURE;
        std::cout << "PASSED" << std::endl << std::endl;

        std::cout << "TEST : stream tuple read/write" << std::endl;
        std::stringstream stream1;
        stream1 << serialize_tuple(tuple1);
        tuple_type tuple3 = read_tuple<tuple_type>(stream1);
        check = tools::is_float_equal(tuple1, tuple2, 1e-6);
        std::cout <<  std::boolalpha << print_tuple(tuple3) << "==" << print_tuple(tuple1) << " = " << check << std::endl;
        std::cout << "PASSED" << std::endl << std::endl;

        std::function<double(int,double,complex_type)> f1 = [](int x, double y, complex_type z){return x+y+real(z);};
        std::cout << "function (" << print_tuple(tuple1) << ") = " << unfold_tuple(f1,tuple1) << std::endl;
        struct struct1 { int operator()(int x,double y,complex_type z){return x+y+real(z);}; };
        std::cout << "structure(" << print_tuple(tuple1) << ") = " << unfold_tuple(struct1(),tuple1) << std::endl;
    }

    {
        auto tuple1 = repeater<double, 4>::get_tuple(3.0);
        std::cout << "tuple:         " << print_tuple(tuple1) << std::endl;
        auto array1 = repeater<double, 4>::get_array(3.0);
        std::cout << "array:         " << print_array(array1) << std::endl;
    }
    return EXIT_SUCCESS;
}
