#ifndef ___FK_GRID_H___
#define ___FK_GRID_H___

#include "Common.h"

namespace GFTools { 

template <typename ValueType, typename ... > struct GridPointTypeExtractor;
template <typename ValueType, typename ... > struct GridPointExtractor;

/** A representation of a one-dimensional grid, which stores an array of the ValueType values. */
template <typename ValueType, class Derived>
class Grid {
public:
    /** A point combines a point of the grid and it's index. */
    struct point {
        ValueType _val;
        size_t _index;
        operator ValueType() const { return _val; }
        explicit inline operator size_t() const { return _index; }
        explicit inline operator int() const { return _index; }
        point(){};
        inline point(ValueType val, size_t index):_val(val),_index(index){};
        inline point(const point& rhs):_val(rhs._val),_index(rhs._index){};
        inline point(point&& rhs) { _val = rhs._val, _index = rhs._index; }
        inline point& operator=(point&& rhs) { _val = rhs._val, _index = rhs._index; return *this;}
        inline point operator=(const point& rhs) { _val = rhs._val, _index = rhs._index; return *this;}
        //inline point(ValueType in):_val(in){_index=2;};
        bool operator==(const point &rhs) const {return (_val == rhs._val) && (_index == rhs._index);}
        friend std::ostream& operator<<(std::ostream& lhs, const point &gr)
            {lhs<<"{"<<gr._val<<"<-["<<gr._index<<"]}"; return lhs;};
    };
protected:
    std::vector<point> _vals;
    Grid(const std::vector<point> & vals);
public:
    typedef typename std::vector<ValueType>::iterator iterator;
    typedef typename std::vector<ValueType>::const_iterator const_iterator;
    typename Grid<ValueType, Derived>::const_iterator begin() const { return _vals.begin(); };
    typename Grid<ValueType, Derived>::const_iterator end() const { return _vals.end(); };
    /** Empty constructor. */
    Grid();
    /** Copy from vector. */
    Grid(const std::vector<ValueType> & vals);
    /** Initialize the values from a given function, that maps the integer values
     * to the ValueType values. 
     */
    Grid(int min, int max, std::function<ValueType (int)> f);
    /** Copy constructor. */
    Grid(const Grid& rhs):_vals(rhs._vals){};
    /** Move constructor. */
    Grid(Grid&& rhs){_vals.swap(rhs._vals);};
    /** Returns a value at given index. */
    point operator[](size_t in) const;
    /** Returns all values. */
    const std::vector<point> & getPoints() const;
    /** Returns values of all points. */
    const std::vector<ValueType> & getPointVals() const;
    /** Returns size of grid. */
    size_t getSize() const;

    /** Shift a point by the given value. */
    template <class ArgType>
        point shift(point in, ArgType shift_arg) const;
    template <class ArgType>
        ValueType shift(ValueType in, ArgType shift_arg) const;

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
    template <class Obj> auto integrate(const Obj &in) const ->decltype(in[_vals[0]]) 
        { return static_cast<const Derived*>(this)->integrate(in); };
    /** Integrate over grid with extra arguments provided. */
    template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...))
        { return static_cast<const Derived*>(this)->integrate(in, Args...); };
    /** Make the object printable. */
    template <typename ValType, class Derived2> friend std::ostream& operator<<(std::ostream& lhs, const Grid<ValType,Derived2> &gr);

    class exWrongIndex : public std::exception { virtual const char* what() const throw(); }; 
};

template <class Grid>
std::ostream& operator<<(std::ostream& lhs, const __num_format<typename Grid::point> &in){lhs << std::setprecision(in._prec) << in._v._val; return lhs;};
template <class Grid>
std::istream& operator>>(std::istream& lhs, __num_format<typename Grid::point> &out)
{__num_format<decltype(std::declval<typename Grid::point>()._v)> d(0.0); lhs >> d; out._v._val = d._v; return lhs;};


/** A tool to generate a function of argtypes of grids. */
template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...GridTypes, typename ...ArgTypes>
struct GridPointTypeExtractor<ValueType, T<GridType1, GridTypes...>, ArgTypes...> : 
GridPointTypeExtractor<ValueType, T<GridTypes...>, ArgTypes...,decltype(GridType1::point::_val)>
{
};

template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...ArgTypes>
struct GridPointTypeExtractor<ValueType, T<GridType1>, ArgTypes...> {
    typedef std::function<ValueType(ArgTypes...,decltype(GridType1::point::_val))> type; 
    typedef std::function<ValueType(ArgTypes...,decltype(GridType1::point::_val))> arg_type; 
    typedef std::tuple<ArgTypes...,decltype(GridType1::point::_val)> arg_tuple_type;
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
};


/* A tool to generate an array of grid sizes from a given tuple of grids. */
template <size_t N>
struct GetGridSizes {
    template <typename... GridType, size_t M>
    static inline void TupleSizeToArray( const std::tuple<GridType...>& in, std::array<size_t, M> &out ) {
        static_assert(N>1,"!");
        std::get<N-1>(out) = std::get<N-1>(in).getSize();
        GetGridSizes<N-1>::TupleSizeToArray( in, out );
    }
};

template <>
template <typename... GridType, size_t M>
inline void GetGridSizes<1>::TupleSizeToArray( const std::tuple<GridType...>& in, std::array<size_t, M> &out ) {
    std::get<0>(out) = std::get<0>(in).getSize();
}


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
inline Grid<ValueType,Derived>::Grid()
{};

template <typename ValueType, class Derived>
inline Grid<ValueType,Derived>::Grid(const std::vector<point> &vals):_vals(vals)
{
};


template <typename ValueType, class Derived>
inline Grid<ValueType,Derived>::Grid(const std::vector<ValueType> &vals):_vals(vals)
{
    for (size_t i=0; i<_vals.size(); ++i) _vals[i]._index = i;
};

template <typename ValueType, class Derived>
inline Grid<ValueType,Derived>::Grid(int min, int max, std::function<ValueType (int)> f)
{
    if (max<min) std::swap(min,max);
    size_t n_points = max-min;
    _vals.resize(n_points); 
    for (int i=0; i<n_points; ++i) _vals[i]= point(f(min+i), i) ; 
}

template <typename ValueType, class Derived>
inline typename Grid<ValueType,Derived>::point Grid<ValueType,Derived>::operator[](size_t index) const
{
    if (index>_vals.size()) throw exWrongIndex();
    return _vals[index];
}

template <typename ValueType, class Derived>
inline const std::vector<typename Grid<ValueType,Derived>::point> & Grid<ValueType,Derived>::getPoints() const
{
    return _vals;
}

template <typename ValueType, class Derived>
inline const std::vector<ValueType>& Grid<ValueType,Derived>::getPointVals() const
{
    return _vals;
}

template <typename ValueType, class Derived>
inline size_t Grid<ValueType,Derived>::getSize() const
{
    return _vals.size();
}

template <typename ValueType, class Derived>
template <class ArgType>
inline typename Grid<ValueType,Derived>::point Grid<ValueType,Derived>::shift(point in, ArgType shift_arg) const
{
    point out;
    out._val = in._val + ValueType(shift_arg);
    auto find_result = this->find(out._val);
    if (std::get<0>(find_result)) { out._index = std::get<1>(find_result); return (*this)[out._index]; }
    else { out._index = this->getSize(); 
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
    return point(_vals[std::get<1>(find_result)],std::get<1>(find_result));
} 

template <typename ValueType, class Derived>
template <class ArgType>
inline ValueType Grid<ValueType,Derived>::shift(ValueType in, ArgType shift_arg) const
{
    return in+ValueType(shift_arg);
}


template <typename ValueType, class Derived>
std::ostream& operator<<(std::ostream& lhs, const Grid<ValueType,Derived> &gr)
{ 
    lhs << "{";
    //lhs << gr._vals;
    std::ostream_iterator<ValueType> out_it (lhs,", ");
    std::transform(gr._vals.begin(),gr._vals.end(), out_it, [](const typename Grid<ValueType,Derived>::point &x){return ValueType(x);});
    lhs << "}";
    return lhs;
}

template <typename ValueType, class Derived>
const char* Grid<ValueType,Derived>::exWrongIndex::what() const throw(){ 
     return "Index out of bounds";
};



} // end :: namespace GFTools
#include "Grid.hpp"

#endif // endin :: ifndef ___FK_GRID_H___
