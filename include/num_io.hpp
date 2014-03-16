#pragma once

#include "defaults.hpp"
#include <fstream>
#include <iomanip>

namespace gftools {

/** num_io is a wrapper around object of type T for various input/output operations, 
    such as streaming, saving to/loading from plaintext files.  */
template <typename T> struct num_io;

/** Read from stream. */
template <typename T> std::ostream& operator<<(std::ostream& lhs, const num_io<T> &in);
/** Write to stream. */
template <typename T> std::istream& operator>>(std::istream& lhs, num_io<T> &out);

template <typename T> 
struct num_io {
    typedef typename std::remove_reference<T>::type type; 

    /// Construct from value
    num_io(type& v):value_(v){};
    /// Delete default constructor
    num_io()=delete;
    /// Delete copy constructor 
    num_io(const num_io&)=delete;

    operator T(){return value_;};
    void savetxt(const std::string& filename) { 
        std::cout << "Saving " << typeid(*this).name() << " to " << filename << std::endl;
        std::ofstream out; out.open(filename.c_str()); out << *this << std::endl; out.close(); 
    }; 
    void loadtxt(const std::string& filename) { 
        std::cout << "Loading " << typeid(*this).name() << " from " << filename << std::endl;
        std::ifstream out; out.open(filename.c_str()); if (out.fail()) throw (std::bad_exception()); out >> *this; out.close(); 
    }; 

    friend std::ostream& operator<< <>(std::ostream& lhs, const num_io<T> &in);
    friend std::istream& operator>> <>(std::istream& lhs, num_io<T> &out);

    /** Underlying value. */
    type& value_;
    /** Output precision. */
    static const int prec_ = 12;
};

template <typename T> 
inline std::ostream& operator<<(std::ostream& lhs, const num_io<T> &in) {lhs << std::setprecision(in.prec_) << in.value_; return lhs;};
template <typename T> 
inline std::istream& operator>>(std::istream& lhs, num_io<T> &out) {lhs >> out.value_; return lhs;};
template <>
inline std::ostream& operator<<(std::ostream& lhs, const num_io<complex_type> &in){lhs << std::setprecision(in.prec_) << real(in.value_) << " " << imag(in.value_); return lhs;};
template <>
inline std::istream& operator>>(std::istream& lhs, num_io<complex_type> &out){real_type re,im; lhs >> re; lhs >> im; out.value_ = re+I*im; return lhs;};


} // end of namespace gftools
