#pragma once

#include "defaults.hpp"
#include <fstream>
#include <iomanip>

namespace gftools {

/** num_io is a wrapper around another object (of type T). 
    num_io does input/output operations.
    such as streaming, saving to/loading from plaintext files.
    Internally it stores references to an instance of type T.
    Typical objects are: numbers, arrays of numbers, complex numbers, grid points*/
template <typename T> class num_io {
public:
    typedef typename std::remove_reference<T>::type type; 

    /// Construct from value
    num_io(type& v):value_(v){}

    /// Access the internally stored value
    type &operator()(){return value_;}
    /// Const access to the internally stored value
    const type &operator()() const{return value_;}
    /// Save the object to a file with a filename passed in as a name
    void savetxt(const std::string& filename) { 
        std::cout << "Saving " << typeid(*this).name() << " to " << filename << std::endl;
        std::ofstream out; out.open(filename.c_str()); out << *this << std::endl; out.close(); 
    }; 
    /// Read the object from a file with a filename passed in as a name
    void loadtxt(const std::string& filename) { 
        std::cout << "Loading " << typeid(*this).name() << " from " << filename << std::endl;
        std::ifstream out; out.open(filename.c_str()); if (out.fail()) throw (std::bad_exception()); out >> *this; out.close(); 
    }; 

    type& value() { return value_; } // TODO : test using it
    type const& value() const { return value_; } // TODO : test using it
    static constexpr int precision() { return prec_; }
    static constexpr double tolerance() { return tol_; }

private:
    /// Underlying value.
    type& value_;
    /// Output precision. TODO: replace c++-11 construct constexpr with something backwards compatible.
    static constexpr int prec_ = std::numeric_limits<double>::max_digits10;
    /// Comparison tolerance
    static constexpr double tol_ = 1e-10;
};

/// A shortcut to make num_io wrapper objects without specifying the template parameter
template <typename T>
num_io<T> make_num_io (T &v){return num_io<T>(v);}

template <typename T> inline std::ostream& operator<<(std::ostream& lhs, const num_io<T> &in) {
  lhs << std::setprecision(in.precision()) << in(); return lhs;
};
template <typename T> 
inline std::istream& operator>>(std::istream& lhs, num_io<T> out) {
  lhs >> out(); return lhs;
}
template <>
inline std::ostream& operator<<(std::ostream& lhs, const num_io<complex_type> &in){lhs << std::setprecision(in.precision()) << real(in()) << " " << imag(in()); return lhs;};
template <>
inline std::istream& operator>>(std::istream& lhs, num_io<complex_type> out){real_type re,im; lhs >> re; lhs >> im; out() = re+I*im; return lhs;};


} // end of namespace gftools
