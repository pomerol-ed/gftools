#pragma once 
#include <iterator>

#include "defaults.hpp"
#include "tools.hpp"

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
    point_base () = delete;
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
        //point():point_base<ValueType>::point_base(){};
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
    /** rvalue assignment. */
    grid_base& operator=(grid_base &&rhs){vals_.swap(rhs.vals_); return (*this);};
    /** Returns a value at given index. */
    point operator[](size_t in) const;
    /** Returns all values. */
    const std::vector<point> & points() const;
    /** Returns values of all points. */
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

    // CFTP forwards
    /** Get a value of an object at the given coordinate, which is defined on a grid. */
    //template <class Obj> 
    //    auto eval(Obj &in, ValueType x) const ->decltype(in[0])
    //    { return static_cast<const Derived*>(this)->eval(in,x); };
    /** Returns a tuple of left closest index, weight, right closest index and weight, which are the closest to input value. */
    /** Integrate over grid. */
    //template <class Obj> auto integrate(const Obj &in) const ->decltype(in[vals_[0]]) 
    //    { return static_cast<const Derived*>(this)->integrate(in); };
    /** Integrate over grid with extra arguments provided. */
    template <class Obj, typename ...OtherArgTypes, typename std::result_of<Obj(point,OtherArgTypes...)>::type> 
        auto integrate(Obj &&in, OtherArgTypes... Args) const -> 
            typename std::remove_reference<typename std::result_of<Obj(point,OtherArgTypes...)>::type>::type
                { return static_cast<const Derived*>(this)->integrate(in, Args...); };
    template <class Obj, decltype (std::declval<Obj>()[std::declval<point>()])>
        auto integrate(Obj &&in) const -> 
            typename std::remove_reference<decltype (std::declval<Obj>()[std::declval<point>()])>::type
                { return static_cast<const Derived*>(this)->integrate(in); };
    /** Make the object printable. */
    template <typename ValType, class Derived2> friend std::ostream& operator<<(std::ostream& lhs, const grid_base<ValType,Derived2> &gr);


    std::vector<point> vals_;
    class ex_wrong_index : public std::exception { virtual const char* what() const throw(); }; 
};

namespace extra { 
/// A small helper struct to decorate a function object with [] method.
template <typename F, typename Grid>
struct function_proxy { 
    F f_;
    const Grid& grid_;
    function_proxy(F f, Grid const& grid):f_(f),grid_(grid){}
    typedef typename std::result_of<F(typename Grid::value_type)>::type value_type;
    value_type operator[](int i) const {return f_((grid_.points()[i]).val_); };
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
    if (index>vals_.size()) throw ex_wrong_index();
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
    for (const auto& x : vals_) out.emplace_back(x.val_);
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
    return (in.index_ < vals_.size() && std::abs(in.val_ - vals_[in.index_].val_) < tolerance);
}

template <typename ValueType, class Derived>
template <class Obj>
inline auto grid_base<ValueType,Derived>::eval(Obj &&in, point x) const ->decltype(in[0]) 
{
    if (check_point(x)) return in[x.index_];
    else { ERROR ("Point not found"); throw ex_wrong_index(); };
}

template <typename ValueType, class Derived>
inline typename grid_base<ValueType,Derived>::point grid_base<ValueType,Derived>::find_nearest(ValueType in) const
{
    static_assert(std::is_same<bool,decltype(std::declval<ValueType>() < std::declval<ValueType>())>::value, 
        "Default find_nearest is written only for less-comparable types");
    auto nearest_iter = std::lower_bound(vals_.begin(), vals_.end(), in, [](ValueType x, ValueType y){return x<y;});
    size_t dist = std::distance(vals_.begin(), nearest_iter);
    if (dist > 0 && std::abs(complex_type(vals_[dist].val_) - complex_type(in)) > std::abs(complex_type(vals_[dist-1].val_) - complex_type(in)) ) dist--;
    return vals_[dist];
} 

template <typename ValueType, class Derived>
inline typename grid_base<ValueType,Derived>::point grid_base<ValueType,Derived>::shift(point in, ValueType shift_arg) const
{
    if (tools::is_float_equal(shift_arg, 0.0)) return in;
    point out(in);
    out.val_ = static_cast<const Derived*>(this)->shift(ValueType(in),shift_arg);
    point p1 = static_cast<const Derived*>(this)->find_nearest(out.val_);
    if (!tools::is_float_equal(p1.val_, out.val_)) { 
        std::cerr << "Couldn't shift point" << std::endl; 
        throw (ex_wrong_index());
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
    size_t index = (in.index_ + shift_arg.index_)%vals_.size();
    #ifndef NDEBUG
    ValueType val = static_cast<const Derived*>(this)->shift(in.val_, shift_arg.val_);
    if (!tools::is_float_equal(val, vals_[index].val_)) throw (ex_wrong_index()); 
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

