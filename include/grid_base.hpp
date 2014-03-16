#ifndef ___GFTOOLS_GRID_H___
#define ___GFTOOLS_GRID_H___

#include <boost/type_traits.hpp> 
#include <boost/utility/enable_if.hpp>
#include "Defaults.hpp"
#include "Tools.hpp"
#include <iterator>

namespace GFTools { 

template <typename ValueType, typename ... > struct GridArgTypeExtractor;
template <typename ValueType, typename ... > struct GridPointExtractor;

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
    friend std::ostream& operator<<(std::ostream& lhs, const point_base &p){lhs<<"{"<<p.val_<<"<-["<<p.index_<<"]}"; return lhs;};

    ValueType val_;
    size_t index_;
}; 

/** A representation of a one-dimensional grid, which stores an array of the ValueType values. */
//template <typename ValueType, class Derived, typename = typename std::enable_if<!std::is_same<ValueType, int>::value>::type>
template <typename ValueType, class Derived>
class Grid {
public:
    /** A point combines a point of the grid and it's index. */
    struct point : point_base<ValueType> { 
        point():point_base<ValueType>::point_base(){};
        point(ValueType val, size_t index):point_base<ValueType>::point_base(val,index){};
        point(const point_base<ValueType> &in):point_base<ValueType>::point_base(in){};
        point(point_base<ValueType> &&in):point_base<ValueType>::point_base(in){};
        };
    typedef ValueType value_type;

    std::vector<point> vals_;

    Grid(const std::vector<point> & vals);
    typedef typename std::vector<ValueType>::iterator iterator;
    typedef typename std::vector<ValueType>::const_iterator const_iterator;
    typename Grid<ValueType, Derived>::const_iterator begin() const { return vals_.begin(); };
    typename Grid<ValueType, Derived>::const_iterator end() const { return vals_.end(); };
    /** Empty constructor. */
    Grid();
    /** Copy from vector. */
    Grid(const std::vector<ValueType> & vals);
    /** Initialize the values from a given function, that maps the integer values
     * to the ValueType values. 
     */
    Grid(int min, int max, std::function<ValueType (int)> f);
    /** Copy constructor. */
    Grid(const Grid& rhs):vals_(rhs.vals_){};
    /** Move constructor. */
    Grid(Grid&& rhs){vals_.swap(rhs.vals_);};
    /** Returns a value at given index. */
    point operator[](size_t in) const;
    /** Returns all values. */
    const std::vector<point> & getPoints() const;
    /** Returns values of all points. */
    const std::vector<ValueType> & getPointVals() const;
    /** Checks if a point is present in a grid. */
    bool checkPoint(point in, RealType tolerance = std::numeric_limits<RealType>::epsilon()) const;
    /** Returns size of grid. */
    size_t getSize() const;

    /** Shift a point by the given value. */
    template <class ArgType>
        point shift(point in, ArgType shift_arg) const;
    template <class ArgType>
        ValueType shift(ValueType in, ArgType shift_arg) const;
    point shift (point in, point shift_arg) const;

    // CFTP forwards
    /** Get a value of an object at the given point, which is defined on a grid. */
    template <class Obj> auto getValue(Obj &in, Grid<ValueType,Derived>::point x) const ->decltype(in[0])
        { return static_cast<const Derived*>(this)->getValue(in,x); };
    /** Get a value of an object at the given coordinate, which is defined on a grid. */
    template <class Obj> auto getValue(Obj &in, ValueType x) const ->decltype(in[0])
        { return static_cast<const Derived*>(this)->getValue(in,x); };
    /** Returns a tuple of left closest index, weight, right closest index and weight, which are the closest to input value. */
    std::tuple <bool, size_t, RealType> find (ValueType in) const 
        { return static_cast<const Derived*>(this)->find(in); };
    /** Returns the closest point to the given value. */
    point findClosest(ValueType in) const;
    /** Integrate over grid. */
    template <class Obj> auto integrate(const Obj &in) const ->decltype(in[vals_[0]]) 
        { return static_cast<const Derived*>(this)->integrate(in); };
    /** Integrate over grid with extra arguments provided. */
    template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(vals_[0],Args...))
        { return static_cast<const Derived*>(this)->integrate(in, Args...); };
    /** Make the object printable. */
    template <typename ValType, class Derived2> friend std::ostream& operator<<(std::ostream& lhs, const Grid<ValType,Derived2> &gr);

    class exWrongIndex : public std::exception { virtual const char* what() const throw(); }; 

};

template <class Grid>
std::ostream& operator<<(std::ostream& lhs, const __num_format<typename Grid::point> &in){lhs << std::setprecision(in._prec) << in._v.val_; return lhs;};
template <class Grid>
std::istream& operator>>(std::istream& lhs, __num_format<typename Grid::point> &out)
{__num_format<decltype(std::declval<typename Grid::point>()._v)> d(0.0); lhs >> d; out._v.val_ = d._v; return lhs;};

/** A tool to generate a function of argtypes of grids. */
template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...GridTypes, typename ...ArgTypes>
struct GridArgTypeExtractor<ValueType, T<GridType1, GridTypes...>, ArgTypes...> : 
GridArgTypeExtractor<ValueType, T<GridTypes...>, ArgTypes...,decltype(GridType1::point::val_)>
{
};

template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...ArgTypes>
struct GridArgTypeExtractor<ValueType, T<GridType1>, ArgTypes...> {
    typedef std::function<ValueType(ArgTypes...,decltype(GridType1::point::val_))> type; 
    typedef std::function<ValueType(ArgTypes...,decltype(GridType1::point::val_))> arg_type; 
    typedef std::tuple<ArgTypes...,decltype(GridType1::point::val_)> arg_tuple_type;
    typedef __caller<ValueType, ArgTypes..., decltype(GridType1::point::val_)> arg_function_wrapper;
};

/** A tool to generate a function of argtypes of grids. */
template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...GridTypes, typename ...ArgTypes>
struct GridPointExtractor<ValueType, T<GridType1, GridTypes...>, ArgTypes...> : 
GridPointExtractor<ValueType, T<GridTypes...>, ArgTypes...,typename GridType1::point>
{
};

template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...ArgTypes>
struct GridPointExtractor<ValueType, T<GridType1>, ArgTypes...> {
    typedef std::function<ValueType(ArgTypes...,typename GridType1::point)> point_type; 
    typedef std::tuple<ArgTypes...,typename GridType1::point> arg_tuple_type;
    typedef __caller<ValueType, ArgTypes..., typename GridType1::point> point_function_wrapper;
};


/* A tool to generate an array of grid sizes from a given tuple of grids. */
template <size_t N> 
struct GetGridSizes {
    template <typename ... GridType>
    static inline std::array<size_t,sizeof...(GridType)> TupleSizeToArray( const std::tuple<GridType...>& in ) {
        static_assert(N>1,"!");
        auto out = GetGridSizes<N-1>::template TupleSizeToArray<GridType...>( in );
        std::get<N-1>(out) = std::get<N-1>(in).getSize();
        return out;
    };
};

template <> 
struct GetGridSizes<1> {
    template <typename... GridType>
    static inline std::array<size_t, sizeof...(GridType)> TupleSizeToArray( const std::tuple<GridType...>& in ) {
        std::array<size_t,sizeof...(GridType)> out;
        std::get<0>(out) = std::get<0>(in).getSize();
        return out;
    }
};

/** A tool to recursiverly integrate over a grid. */
template <typename GridType, class Obj> struct RecursiveGridIntegrator;

/* Integrate a function over a grid. */
template <typename GridType, typename ValueType, typename ArgType1> 
struct RecursiveGridIntegrator<GridType, ValueType(ArgType1)>
{
    inline static ValueType integrate(const GridType& grid, const std::function<ValueType(ArgType1)>& in){ 
        //DEBUG("Using this method, N=1 with "<< sizeof...(OtherArgs) << " other args" );
        return grid.integrate(in);
    }
};

/** An alias for a std::function template type. */
template <typename GridType, typename ValueType, typename ArgType1> 
struct RecursiveGridIntegrator<GridType, std::function<ValueType(ArgType1)>>:RecursiveGridIntegrator<GridType, ValueType(ArgType1)>{};

/** Multi-dimensional integration over the same grid. */
template <typename GridType, typename ValueType, typename ArgType1, typename ...ArgTypes> 
struct RecursiveGridIntegrator<GridType, ValueType(ArgType1, ArgTypes...)>
{
    typedef std::function<ValueType(ArgTypes...)> type;
    inline static ValueType integrate(
        const GridType& grid, 
        const std::function<ValueType(ArgType1, ArgTypes...)>& in) 
     {
        //DEBUG("Using this method, N!=1");
        auto f1 = [&](const ArgType1& arg1)->ValueType { 
            //DEBUG(arg1);
            std::function<ValueType(ArgTypes...)> f0 = [&](const ArgTypes&... args){return in(arg1,args...);};
            return RecursiveGridIntegrator<GridType, ValueType(ArgTypes...)>::integrate(grid, f0);
            };
        //DEBUG(f1(0.0));
        return grid.integrate(f1);
    } 
};

/** An alias for a std::function template type. */
template <typename GridType, typename ValueType, typename ArgType1, typename ...ArgTypes> 
struct RecursiveGridIntegrator<GridType, std::function<ValueType(ArgType1, ArgTypes...)>>:
    RecursiveGridIntegrator<GridType, ValueType(ArgType1, ArgTypes...)> {};

//
// Grid implementation
//

template <typename ValueType, class Derived>
Grid<ValueType,Derived>::Grid()
{};

template <typename ValueType, class Derived>
Grid<ValueType,Derived>::Grid(const std::vector<point> &vals):vals_(vals)
{
};


template <typename ValueType, class Derived>
Grid<ValueType,Derived>::Grid(const std::vector<ValueType> &vals)
{
    vals_.reserve(vals.size());
    for (size_t i=0; i<vals_.size(); ++i) { vals_.emplace_back(point(i, vals[i])); };
};

template <typename ValueType, class Derived>
Grid<ValueType,Derived>::Grid(int min, int max, std::function<ValueType (int)> f)
{
    if (max<min) std::swap(min,max);
    size_t n_points = max-min;
    vals_.reserve(n_points);
    for (int i=0; i<n_points; ++i) vals_.emplace_back(f(min+i), i) ; 
}

template <typename ValueType, class Derived>
inline typename Grid<ValueType,Derived>::point Grid<ValueType,Derived>::operator[](size_t index) const
{
    if (index>vals_.size()) throw exWrongIndex();
    return vals_[index];
}

template <typename ValueType, class Derived>
inline const std::vector<typename Grid<ValueType,Derived>::point> & Grid<ValueType,Derived>::getPoints() const
{
    return vals_;
}

template <typename ValueType, class Derived>
inline const std::vector<ValueType>& Grid<ValueType,Derived>::getPointVals() const
{
    return vals_;
}

template <typename ValueType, class Derived>
inline size_t Grid<ValueType,Derived>::getSize() const
{
    return vals_.size();
}

template <typename ValueType, class Derived>
inline bool Grid<ValueType,Derived>::checkPoint(point in, RealType tolerance) const
{
    return (in.index_ < vals_.size() && std::abs(in.val_ - vals_[in.index_].val_) < tolerance);
}


template <typename ValueType, class Derived>
template <class ArgType>
inline typename Grid<ValueType,Derived>::point Grid<ValueType,Derived>::shift(point in, ArgType shift_arg) const
{
    point out;
    if (std::abs(ValueType(shift_arg))<std::numeric_limits<RealType>::epsilon()) return in;
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
inline typename Grid<ValueType,Derived>::point Grid<ValueType,Derived>::findClosest(ValueType in) const
{
    auto find_result = find(in);
    if (!std::get<0>(find_result)) { ERROR("Couldn't find the closest point"); throw (exWrongIndex()); };
    return point(vals_[std::get<1>(find_result)],std::get<1>(find_result));
} 

template <typename ValueType, class Derived>
template <class ArgType>
inline ValueType Grid<ValueType,Derived>::shift(ValueType in, ArgType shift_arg) const
{
    return in+ValueType(shift_arg);
}

template <typename ValueType, class Derived>
inline typename Grid<ValueType,Derived>::point Grid<ValueType,Derived>::shift(point in, point shift_arg) const
{
    size_t index = (in.index_ + shift_arg.index_)%vals_.size();
    #ifndef NDEBUG
    ValueType val = this->shift(in.val_, shift_arg.val_);
    if (std::abs(val - vals_[index].val_)>1e-3) throw (exWrongIndex()); 
    #endif
    return vals_[index];

}

template <typename ValueType, class Derived>
std::ostream& operator<<(std::ostream& lhs, const Grid<ValueType,Derived> &gr)
{ 
    lhs << "{";
    //lhs << gr.vals_;
    std::ostream_iterator<ValueType> out_it (lhs,", ");
    std::transform(gr.vals_.begin(),gr.vals_.end(), out_it, [](const typename Grid<ValueType,Derived>::point &x){return ValueType(x);});
    lhs << "}";
    return lhs;
}

template <typename ValueType, class Derived>
const char* Grid<ValueType,Derived>::exWrongIndex::what() const throw(){ 
     return "Index out of bounds";
};

} // end :: namespace GFTools

#endif // endin :: ifndef ___GFTOOLS_GRID_H___
