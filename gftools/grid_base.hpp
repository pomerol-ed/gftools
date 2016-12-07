#pragma once 

#include <boost/operators.hpp>

#include <iterator>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

#include "defaults.hpp"
#include "exceptions.hpp"
#include "num_io.hpp"
#include "almost_equal.hpp"


namespace gftools { 

/** This class describes a point on a grid. It has an index (integer) and a value (pretty much anything) , and it 
maps between the two. The index is used for fast access, and the value is used for other things like interpolation/physics/math.*/
template <typename ValueType> class point_base : 
    ///dependence on this provides additional comparison operators
    boost::less_than_comparable<point_base<ValueType> >,
    ///dependence on this provides additional comparison operators
    boost::equality_comparable<point_base<ValueType> >
 {
public:
    //there is a small wrapper for integers down below. Here we exclude grids of ints to avoid confusion in the cast operators. 
    static_assert(!std::is_same<ValueType,int>::value, "Can't create a grid of ints");

    typedef ValueType value_type;

    ///cast operator to value type
    operator ValueType() const { return val_; }
    ///cast operator to index type
    explicit operator size_t() const { return index_; }
    ///another cast operator to index type
    explicit operator int() const { return index_; }

    ///constructor with a pair of values and indices
    point_base(ValueType val, size_t index):val_(val),index_(index){}
    
    bool operator==(const point_base &rhs) const {return index_ == rhs.index_;}
    bool operator<(const point_base &rhs) const {return this->index_ < rhs.index_;}

    ValueType value() const { return val_; }
    size_t index() const { return index_; }

protected:
    ///grid point (in physical units)
    ValueType val_;
    ///grid point index
    size_t index_;
}; 
    
template<typename T> std::ostream& operator<<(std::ostream& lhs, const point_base<T> &p){
  lhs<<"{"<<p.value()<<"<-["<<p.index()<<"]}"; return lhs;
}

/** A one-dimensional grid, which stores an array of values. Typical examples: grid of real frequencies. Grid of k-points. Grid of Matsubara frequencies. Grid of imaginary times.*/
template <typename ValueType, class Derived>
class grid_base : public boost::equality_comparable<grid_base<ValueType, Derived> > {
public:
    typedef point_base<ValueType> point;
    typedef ValueType value_type;

    /** constructor a grid out of thin air. */
    grid_base();
    /** construct a grid given a vector of points. */
    grid_base(const std::vector<point> & vals);
    /** construct a grid from a vector of values (but not associated indices that would be stored in points). */
    grid_base(const std::vector<ValueType> & vals);
    /** Initialize the values from an external function that maps the integer values to the ValueType values. */
    grid_base(int min, int max, std::function<ValueType (int)> f);
    
    /** Returns a value at given index. */
    point operator[](size_t in) const;
    /** Returns all values. */
    const std::vector<point> & points() const;
    /** Returns values of all points. (computed on the fly, slow)*/
    std::vector<ValueType> values() const;
    /** Checks if a point is present in a grid. */
    bool check_point(point in, real_type tolerance = std::numeric_limits<real_type>::epsilon()) const;
    /** Returns size of grid. */
    size_t size() const;

    /** Returns the closest point to the given value. */
    point find_nearest(ValueType in) const;
    /** Get a value of an object at the given point, which is defined on a grid. */
    template <class Obj> auto eval(Obj &&in, point x) const ->decltype(in[0]);

    /** Shift a point by the given value. */
    point shift(point in, ValueType shift_arg) const;
    ValueType shift(ValueType in, ValueType shift_arg) const;
    point shift (point in, point shift_arg) const;

    // CRTP forwards: TODO: these will need to be cleaned for C++ and documented.
    /** Get a value of an object at the given coordinate, which is defined on a grid. */
    template <class Obj>
        auto eval(Obj &in, ValueType x) const ->decltype(in[0])
        { return static_cast<const Derived*>(this)->eval(in,x); };
    /// Make the object printable.
    template <typename ValType, class Derived2> friend std::ostream& operator<<(std::ostream& lhs, const grid_base<ValType,Derived2> &gr);
    /// Compare 2 grids
    bool operator==(const grid_base &rhs) const;

    class ex_wrong_index : public gftools_exception { public: 
        ex_wrong_index(int i, int l):index_(i),l_(l){}; 
        virtual const char* what() const throw(){return std::string("grid_base : index " + std::to_string(index_) + "is out of bounds >" + std::to_string(l_) + ".").c_str();} 
        int index_; int l_; 
        };

    class ex_not_found : public gftools_exception { public: 
        ex_not_found(value_type x, grid_base const &y):x_(x), y_(y){}; 
        virtual const char* what() const throw(){value_type b1(y_[0]), b2(y_[y_.size()-1]); return std::string("grid_base : " + make_num_io(x_).to_string() + " is not found in \
the grid [" + make_num_io(b1).to_string() + "; " + make_num_io(b2).to_string() + "]." ).c_str();} 
        value_type x_; 
        grid_base const& y_;
        };


protected:
    std::vector<point> vals_;
};

namespace extra { 
/// A small helper struct to decorate a function object with [] method.
template <typename F, typename Grid>
struct function_proxy { 
    F f_;
    const Grid& grid_;
    function_proxy(F f, Grid const& grid):f_(f),grid_(grid){}
    typedef typename std::result_of<F(typename Grid::value_type)>::type value_type;
    value_type operator[](int i) const {return f_((grid_.points()[i]).value()); }
};
}



//
// grid_base implementation
//

template <typename ValueType, class Derived>
grid_base<ValueType,Derived>::grid_base()
{};

template <typename ValueType, class Derived>
grid_base<ValueType,Derived>::grid_base(const std::vector<point> &vals):vals_(vals)
{
};


template <typename ValueType, class Derived>
grid_base<ValueType,Derived>::grid_base(const std::vector<ValueType> &vals)
{
    vals_.reserve(vals.size());
    for (size_t i=0; i<vals.size(); ++i) { vals_.emplace_back(point(vals[i], i)); };
};

template <typename ValueType, class Derived>
grid_base<ValueType,Derived>::grid_base(int min, int max, std::function<ValueType (int)> f)
{
    if (max<min) std::swap(min,max);
    size_t n_points = max-min;
    vals_.reserve(n_points);
    for (int i=0; i<n_points; ++i) vals_.emplace_back(f(min+i), i) ; 
}

template <typename ValueType, class Derived>
inline typename grid_base<ValueType,Derived>::point grid_base<ValueType,Derived>::operator[](size_t index) const
{
    #ifndef NDEBUG
    if (index>vals_.size()) throw ex_wrong_index(index,this->size()); 
    #endif
    return vals_[index];
}

template <typename ValueType, class Derived>
inline const std::vector<typename grid_base<ValueType,Derived>::point> & grid_base<ValueType,Derived>::points() const
{
    return vals_;
}

template <typename ValueType, class Derived>
inline std::vector<ValueType> grid_base<ValueType,Derived>::values() const
{
    std::vector<ValueType> out;
    out.reserve(vals_.size());
    for (const auto& x : vals_) out.emplace_back(x.value());
    return out;
}

template <typename ValueType, class Derived>
inline size_t grid_base<ValueType,Derived>::size() const
{
    return vals_.size();
}

template <typename ValueType, class Derived>
inline bool grid_base<ValueType,Derived>::check_point(point in, real_type tolerance) const
{
    return (in.index() < vals_.size() && std::abs(in.value() - vals_[in.index()].value()) < tolerance);
}

template <typename ValueType, class Derived>
template <class Obj>
inline auto grid_base<ValueType,Derived>::eval(Obj &&in, point x) const ->decltype(in[0]) 
{
    if (check_point(x)) return in[x.index()];
    else throw ex_wrong_index(x.index(), this->size());
}

template <typename ValueType, class Derived>
inline typename grid_base<ValueType,Derived>::point grid_base<ValueType,Derived>::find_nearest(ValueType in) const
{
    static_assert(std::is_same<bool,decltype(std::declval<ValueType>() < std::declval<ValueType>())>::value, 
        "Default find_nearest is written only for less-comparable types");
    auto nearest_iter = std::lower_bound(vals_.begin(), vals_.end(), in, [](ValueType x, ValueType y){return x<y;});
    size_t dist = std::distance(vals_.begin(), nearest_iter);
    if (dist > 0 && std::abs(complex_type(vals_[dist].value()) - complex_type(in)) > std::abs(complex_type(vals_[dist-1].value()) - complex_type(in)) ) dist--;
    return vals_[dist];
} 

template <typename ValueType, class Derived>
inline typename grid_base<ValueType,Derived>::point grid_base<ValueType,Derived>::shift(point in, ValueType shift_arg) const
{
    if (almost_equal(shift_arg, 0.0)) return in;
    ValueType out(in.value());
    out = static_cast<const Derived*>(this)->shift(ValueType(in),shift_arg);
    point p1 = static_cast<const Derived*>(this)->find_nearest(out);
    if (!almost_equal(p1.value(), out, std::abs(p1.value() - ((p1.index()!=0)?vals_[p1.index() - 1]:vals_[p1.index()+1]).value())/10.)) { 
#ifndef NDEBUG
      ERROR("Couldn't shift point" <<  in << " by " << shift_arg << " got " << out);
#endif
      throw gftools::ex_generic("grid_base::shift : Couldn't shift point");
    }
    else return p1;
}

template <typename ValueType, class Derived>
inline ValueType grid_base<ValueType,Derived>::shift(ValueType in, ValueType shift_arg) const
{
    return in+shift_arg;
}

template <typename ValueType, class Derived>
inline typename grid_base<ValueType,Derived>::point grid_base<ValueType,Derived>::shift(point in, point shift_arg) const
{
    size_t index = (in.index() + shift_arg.index())%vals_.size();
    #ifndef NDEBUG
    ValueType val = static_cast<const Derived*>(this)->shift(in.value(), shift_arg.value());
    if (!almost_equal(val, vals_[index].value())) throw gftools::ex_generic("grid_base::shift : almost equal failed."); 
    #endif
    return vals_[index];

}

template <typename ValueType, class Derived>
inline bool grid_base<ValueType,Derived>::operator==(const grid_base &rhs) const 
{ 
    bool out = (this->size() == rhs.size()); 
    for (int i=0; i<vals_.size() && out; i++) { 
        out = out && almost_equal(vals_[i].value(), rhs.vals_[i].value(), num_io<double>::tolerance()); 
    }
    return out;
}


template <typename ValueType, class Derived>
std::ostream& operator<<(std::ostream& lhs, const grid_base<ValueType,Derived> &gr)
{ 
    lhs << "{";
    //lhs << gr.vals_;
    std::ostream_iterator<ValueType> out_it (lhs,", ");
    std::transform(gr.vals_.begin(),gr.vals_.end(), out_it, [](const typename grid_base<ValueType,Derived>::point &x){return ValueType(x);});
    lhs << "}";
    return lhs;
}

} // end :: namespace gftools
