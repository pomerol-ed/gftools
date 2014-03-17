#ifndef ___GFTOOLS_GRID_H___
#define ___GFTOOLS_GRID_H___

#include "defaults.hpp"
#include "num_io.hpp"
#include "tools.hpp"
//#include "Tools.hpp"
//#include <iterator>

namespace gftools { 

template <typename ValueType>
struct point_base {
    static_assert(!std::is_same<ValueType,int>::value, "Can't create a grid of ints");

    typedef ValueType value_type;

    operator ValueType() const { return val_; }
    explicit operator size_t() const { return index_; }
    explicit operator int() const { return index_; }

    point_base(ValueType val, size_t index):val_(val),index_(index){};
    point_base(const point_base& rhs):val_(rhs.val_),index_(rhs.index_){};
    point_base(point_base&& rhs):val_(rhs.val_),index_(rhs.index_) {};
    point_base& operator=(point_base&& rhs) { val_ = rhs.val_, index_ = rhs.index_; return *this;}
    point_base operator=(const point_base& rhs) { val_ = rhs.val_, index_ = rhs.index_; return *this;}
    point_base (){};
    bool operator==(const point_base &rhs) const {return (val_ == rhs.val_) && (index_ == rhs.index_);}
    bool operator!=(const point_base &rhs) const {return !(*this == rhs); }
    friend std::ostream& operator<<(std::ostream& lhs, const point_base &p){lhs<<"{"<<p.val_<<"<-["<<p.index_<<"]}"; return lhs;};

    ValueType val_;
    size_t index_;
}; 

/** A representation of a one-dimensional grid, which stores an array of the ValueType values. */
//template <typename ValueType, class Derived, typename = typename std::enable_if<!std::is_same<ValueType, int>::value>::type>
template <typename ValueType, class Derived>
class grid_base {
public:
    /** A point combines a point of the grid and it's index. */
    struct point : point_base<ValueType> { 
        point():point_base<ValueType>::point_base(){};
        point(ValueType val, size_t index):point_base<ValueType>::point_base(val,index){};
        point(const point_base<ValueType> &in):point_base<ValueType>::point_base(in){};
        point(point_base<ValueType> &&in):point_base<ValueType>::point_base(in){};
        };

    typedef ValueType value_type;


    grid_base(const std::vector<point> & vals);

    /** Empty constructor. */
    grid_base();
    /** Copy from vector. */
    grid_base(const std::vector<ValueType> & vals);
    /** Initialize the values from a given function, that maps the integer values
     * to the ValueType values. 
     */
    grid_base(int min, int max, std::function<ValueType (int)> f);
    /** Copy constructor. */
    grid_base(const grid_base& rhs):vals_(rhs.vals_){};
    /** Move constructor. */
    grid_base(grid_base&& rhs){vals_.swap(rhs.vals_);};
    /** Returns a value at given index. */
    point operator[](size_t in) const;
    /** Returns all values. */
    const std::vector<point> & get_points() const;
    /** Returns values of all points. */
    const std::vector<ValueType> & get_values() const;
    /** Checks if a point is present in a grid. */
    bool check_point(point in, real_type tolerance = std::numeric_limits<real_type>::epsilon()) const;
    /** Returns size of grid. */
    size_t size() const;

    /** Get a value of an object at the given point, which is defined on a grid. */
    template <class Obj> auto evaluate(Obj &in, point x) const ->decltype(in[0]);

    /** Shift a point by the given value. */
    template <class ArgType>
        point shift(point in, ArgType shift_arg) const;
    template <class ArgType>
        ValueType shift(ValueType in, ArgType shift_arg) const;
    point shift (point in, point shift_arg) const;

    // CFTP forwards
    /** Get a value of an object at the given coordinate, which is defined on a grid. */
    //template <class Obj> 
    //    auto evaluate(Obj &in, ValueType x) const ->decltype(in[0])
    //    { return static_cast<const Derived*>(this)->evaluate(in,x); };
    /** Returns a tuple of left closest index, weight, right closest index and weight, which are the closest to input value. */
    std::tuple <bool, size_t, real_type> find (ValueType in) const 
        { return static_cast<const Derived*>(this)->find(in); };
    /** Returns the closest point to the given value. */
    point find_closest(ValueType in) const;
    /** Integrate over grid. */
    //template <class Obj> auto integrate(const Obj &in) const ->decltype(in[vals_[0]]) 
    //    { return static_cast<const Derived*>(this)->integrate(in); };
    /** Integrate over grid with extra arguments provided. */
    //template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(vals_[0],Args...))
    //    { return static_cast<const Derived*>(this)->integrate(in, Args...); };
    /** Make the object printable. */
    template <typename ValType, class Derived2> friend std::ostream& operator<<(std::ostream& lhs, const grid_base<ValType,Derived2> &gr);


    std::vector<point> vals_;
    class ex_wrong_index : public std::exception { virtual const char* what() const throw(); }; 
};

template <class Grid>
std::ostream& operator<<(std::ostream& lhs, const num_io<typename Grid::point> &in){lhs << std::setprecision(in.prec_) << in.value_.val_; return lhs;};
template <class Grid>
std::istream& operator>>(std::istream& lhs, num_io<typename Grid::point> &out)
{num_io<decltype(std::declval<typename Grid::point>().value_)> d(0.0); lhs >> d; out.value_.val_ = d.value_; return lhs;};

//
// Grid implementation
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
    if (index>vals_.size()) throw ex_wrong_index();
    return vals_[index];
}

template <typename ValueType, class Derived>
inline const std::vector<typename grid_base<ValueType,Derived>::point> & grid_base<ValueType,Derived>::get_points() const
{
    return vals_;
}

template <typename ValueType, class Derived>
inline const std::vector<ValueType>& grid_base<ValueType,Derived>::get_values() const
{
    return vals_;
}

template <typename ValueType, class Derived>
inline size_t grid_base<ValueType,Derived>::size() const
{
    return vals_.size();
}

template <typename ValueType, class Derived>
inline bool grid_base<ValueType,Derived>::check_point(point in, real_type tolerance) const
{
    return (in.index_ < vals_.size() && std::abs(in.val_ - vals_[in.index_].val_) < tolerance);
}

template <typename ValueType, class Derived>
template <class Obj>
inline auto grid_base<ValueType,Derived>::evaluate(Obj &in, point x) const ->decltype(in[0]) 
{
    if (check_point(x)) return in[x.index_];
    else { ERROR ("Point not found"); throw ex_wrong_index(); };
}


template <typename ValueType, class Derived>
template <class ArgType>
inline typename grid_base<ValueType,Derived>::point grid_base<ValueType,Derived>::shift(point in, ArgType shift_arg) const
{
    point out;
    if (std::abs(ValueType(shift_arg))<std::numeric_limits<real_type>::epsilon()) return in;
    out.val_ = in.val_ + ValueType(shift_arg);
    auto find_result = this->find(out.val_);
    if (std::get<0>(find_result)) { out.index_ = std::get<1>(find_result); return (*this)[out.index_]; }
    else { out.index_ = this->getSize(); 
           #ifndef NDEBUG
           ERROR("Returning point with an invalid index after shift.");
           #endif
           return out; };
}

template <typename ValueType, class Derived>
inline typename grid_base<ValueType,Derived>::point grid_base<ValueType,Derived>::find_closest(ValueType in) const
{
    auto find_result = this->find(in);
    if (!std::get<0>(find_result)) { ERROR("Couldn't find the closest point"); throw (ex_wrong_index()); };
    return point(vals_[std::get<1>(find_result)],std::get<1>(find_result));
} 

template <typename ValueType, class Derived>
template <class ArgType>
inline ValueType grid_base<ValueType,Derived>::shift(ValueType in, ArgType shift_arg) const
{
    return in+ValueType(shift_arg);
}

template <typename ValueType, class Derived>
inline typename grid_base<ValueType,Derived>::point grid_base<ValueType,Derived>::shift(point in, point shift_arg) const
{
    size_t index = (in.index_ + shift_arg.index_)%vals_.size();
    #ifndef NDEBUG
    ValueType val = this->shift(in.val_, shift_arg.val_);
    if (std::abs(val - vals_[index].val_)>1e-3) throw (ex_wrong_index()); 
    #endif
    return vals_[index];

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

template <typename ValueType, class Derived>
const char* grid_base<ValueType,Derived>::ex_wrong_index::what() const throw(){ 
     return "Index out of bounds";
};

} // end :: namespace gftools

//#include "grid_tools.hpp"

#endif // endin :: ifndef ___GFTOOLS_GRID_H___
