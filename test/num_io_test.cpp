#include "defaults.hpp"
#include "num_io.hpp"
#include "tools.hpp"

using namespace gftools;

int main()
{
    { // IO of double
        std::cout << "TEST: IO for doubles. " << std::endl;
        typedef double test_type;
        test_type a = 2.147925020934;
        std::cout << "(double) a = " << a << std::endl;
        num_io<test_type> a_io(a);
        a_io.savetxt("a.dat");
        
        test_type b = 1.0;
        num_io<test_type> b_io(b);
        b_io.loadtxt("a.dat");
        std::cout << "(double) b == a = " << std::boolalpha <<  tools::is_float_equal(a,b,a_io.precision()) << std::endl;
        if (!tools::is_float_equal(a,b,a_io.precision())) return EXIT_FAILURE;

        test_type c = b_io();
        std::cout << "(double) c == a = " << std::boolalpha <<  tools::is_float_equal(a,c,a_io.precision()) << std::endl;
        if (!tools::is_float_equal(a,c,a_io.precision())) return EXIT_FAILURE;
        std::cout << "PASSED : " << a << std::endl << std::endl;
    };

    { // IO of complex 
        typedef complex_type test_type;
        std::cout << "TEST: IO for complex numbers. " << std::endl;
        test_type a = 2391.101203+I*8.201007;
        std::cout << "a = " << a << std::endl;
        num_io<decltype(a)> a_io(a);
        a_io.savetxt("a.dat");
        
        test_type b = 1.0;
        num_io<decltype(b)> b_io(b);
        b_io.loadtxt("a.dat");
        std::cout << "b == a = " << std::boolalpha <<  tools::is_float_equal(a,b,a_io.precision()) << std::endl;
        if (!tools::is_float_equal(a,b,a_io.precision())) return EXIT_FAILURE;

        test_type c = b_io();
        std::cout << "c == a = " << std::boolalpha <<  tools::is_float_equal(a,c,a_io.precision()) << std::endl;
        if (!tools::is_float_equal(a,c,a_io.precision())) return EXIT_FAILURE;
        std::cout << "PASSED : " << a << std::endl << std::endl;
    };

    { // IO of int 
        std::cout << "TEST: IO for int. " << std::endl;
        typedef int test_type;
        test_type a = 9179;
        std::cout << "a = " << a << std::endl;
        num_io<decltype(a)> a_io(a);
        a_io.savetxt("a.dat");
        
        test_type b;
        num_io<decltype(b)> b_io(b);
        b_io.loadtxt("a.dat");
        std::cout << "b == a = " << std::boolalpha <<  tools::is_float_equal(a,b,a_io.precision()) << std::endl;
        if (!tools::is_float_equal(a,b,a_io.precision())) return EXIT_FAILURE;

        test_type c = b_io();
        std::cout << "c == a = " << std::boolalpha <<  tools::is_float_equal(a,c,a_io.precision()) << std::endl;
        if (!tools::is_float_equal(a,c,a_io.precision())) return EXIT_FAILURE;
        std::cout << "PASSED : " << a << std::endl << std::endl;
    };



    return EXIT_SUCCESS;
}
