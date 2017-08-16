#pragma once

#define BOOST_RESULT_OF_USE_DECLTYPE

#include<complex>
#include<iostream>
#include<utility>
#include<typeinfo>
#include<boost/core/demangle.hpp>


namespace gftools {

// DEBUG messages with custom verbosity.
// Adapted from http://efesx.com/2010/08/31/overloading-macros/
#ifndef NDEBUG
#define MSG_PREFIX            __FILE__ << ":" << __LINE__ << ": "
#define DEBUG3(MSG,VERBOSITY,VERB_LEVEL)            if (VERBOSITY >= VERB_LEVEL) std::cerr << MSG_PREFIX << MSG << std::endl;
#else
#define DEBUG3(MSG,VERBOSITY,VERB_LEVEL)            ;
#endif

#define DEBUG1(MSG) DEBUG3(MSG,3,3)
#define DEBUG2(MSG,VERBOSITY) DEBUG3(MSG,VERBOSITY,3)

#define VA_NUM_ARGS_IMPL(_1,_2,_3,_4,_5,N,...) N
#define VA_NUM_ARGS(...) VA_NUM_ARGS_IMPL(__VA_ARGS__, 5,4,3,2,1)
#define macro_dispatcher(func, ...) \
            macro_dispatcher_(func, VA_NUM_ARGS(__VA_ARGS__))
#define macro_dispatcher_(func, nargs) \
            macro_dispatcher__(func, nargs)
#define macro_dispatcher__(func, nargs) \
            func ## nargs

#define DEBUG(...) macro_dispatcher(DEBUG, __VA_ARGS__)(__VA_ARGS__)

// some more output presets
#ifndef INFO
#define INFO(MSG)             std::cout << MSG << std::endl;
#define INFO_NONEWLINE(MSG)   std::cout << MSG << std::flush;
#define INFO2(MSG)            std::cout << "\t"     << MSG << std::endl;
#define INFO3(MSG)            std::cout << "\t\t"   << MSG << std::endl;
#define INFO4(MSG)            std::cout << "\t\t\t" << MSG << std::endl;
#endif

#ifndef ERROR
#define MSG_PREFIX            __FILE__ << ":" << __LINE__ << ": "
#define ERROR(MSG)            std::cerr << MSG_PREFIX << MSG << std::endl;
#endif

/** Extract demangled name from a type_info object */
std::string demangled_name(std::type_info const& info) {
    return boost::core::demangle(info.name());
}

typedef double real_type;
typedef std::complex<real_type> complex_type;

/** A short name for imaginary unit. TODO: remove this in dependent codes.*/
static const complex_type I = complex_type(0.0,1.0);    // 'static' to prevent linking problems
/** A short name for pi. TODO: remove this in dependent codes.*/
const real_type PI = M_PI;

} // end namespace gftools

