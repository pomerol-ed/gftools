#pragma once

#include "grid_base.hpp"
#include "tuple_tools.hpp"

namespace gftools {

namespace tools { 

/** A tool to extract ArgTypes(types of numbers) from a tuple of grids. */
template <typename ValueType, typename ... > struct GridArgTypeExtractor;
/** A tool to extract point types from a tuple of grids. */
template <typename ValueType, typename ... > struct GridPointExtractor;

/** A tuple of grids type traits. */
template <typename> struct grid_tuple_traits;

template <typename ... GridTypes>
struct grid_tuple_traits<std::tuple<GridTypes...> >
{
    constexpr static size_t N = sizeof...(GridTypes);
    typedef std::tuple<GridTypes...> grid_tuple_type;

    /** A typedef for a tuple of grid points. */
    typedef typename GridPointExtractor<std::true_type, std::tuple<GridTypes...> >::arg_tuple point_tuple;
    /** A typedef for a tuple of grid point values. */
    typedef typename GridArgTypeExtractor<std::true_type, std::tuple<GridTypes...> >::arg_tuple arg_tuple;
    /** A typedef for a set of indices. */
    typedef std::array<size_t, N> indices;

    template <int M = N-1, typename std::enable_if<M >= 1, bool>::type = 0> 
        static point_tuple get_points(indices in, const grid_tuple_type& grids);
    template <int M = 0,   typename std::enable_if<M == 0, bool>::type = 0> 
        static point_tuple get_points(indices in, const grid_tuple_type& grids);
    static point_tuple get_points(indices in, const grid_tuple_type& grids);

    template <int M = N-1, typename std::enable_if<M >= 1, bool>::type = 0> 
        static arg_tuple get_args(indices in, const grid_tuple_type& grids);
    template <int M = 0,   typename std::enable_if<M == 0, bool>::type = 0> 
        static arg_tuple get_args(indices in, const grid_tuple_type& grids);
    static arg_tuple get_args(indices in, const grid_tuple_type& grids);

    template <int M = N-1, typename std::enable_if<M >= 1, bool>::type = 0> 
        static indices get_indices(point_tuple in, const grid_tuple_type& grids);
    template <int M = 0,   typename std::enable_if<M == 0, bool>::type = 0> 
        static indices get_indices(point_tuple in, const grid_tuple_type& grids);
    static indices get_indices(point_tuple in, const grid_tuple_type& grids);
};


template <typename ...GridTypes>
inline typename grid_tuple_traits<std::tuple<GridTypes...>>::arg_tuple 
    grid_tuple_traits<std::tuple<GridTypes...>>::get_args(indices in, const grid_tuple_type& grids)
{
    return get_args<N-1>(in,grids);
}

template <typename ...GridTypes>
template <int M, typename std::enable_if<M ==0, bool>::type>
inline typename grid_tuple_traits<std::tuple<GridTypes...>>::arg_tuple 
    grid_tuple_traits<std::tuple<GridTypes...>>::get_args(indices in, const grid_tuple_type& grids)
{
    arg_tuple out;
    auto t1 = std::get<N-1>(grids)[in[N-1]].val_;
    std::get<N-1>(out)=t1;
    return out;
}

template <typename ...GridTypes>
template <int M, typename std::enable_if<M >= 1, bool>::type >
inline typename grid_tuple_traits<std::tuple<GridTypes...>>::arg_tuple 
    grid_tuple_traits<std::tuple<GridTypes...>>::get_args(indices in, const grid_tuple_type& grids)
{
    auto out = get_args<M-1>(in,grids);
    auto t1 = std::get<N-1-M>(grids)[in[N-1-M]];
    std::get<N-1-M>(out) = t1.val_;
    return out;
}

template <typename ...GridTypes>
inline typename grid_tuple_traits<std::tuple<GridTypes...>>::point_tuple 
    grid_tuple_traits<std::tuple<GridTypes...>>::get_points(indices in, const grid_tuple_type& grids)
{
    return get_points<N-1>(in,grids);
}

template <typename ...GridTypes>
template <int M, typename std::enable_if<M ==0, bool>::type>
inline typename grid_tuple_traits<std::tuple<GridTypes...>>::point_tuple 
    grid_tuple_traits<std::tuple<GridTypes...>>::get_points(indices in, const grid_tuple_type& grids)
{
    point_tuple out;
    auto t1 = std::get<N-1>(grids)[in[N-1]];
    std::get<N-1>(out)=t1;
    return out;
}

template <typename ...GridTypes>
template <int M, typename std::enable_if<M >= 1, bool>::type >
inline typename grid_tuple_traits<std::tuple<GridTypes...>>::point_tuple 
    grid_tuple_traits<std::tuple<GridTypes...>>::get_points(indices in, const grid_tuple_type& grids)
{
    auto out = get_points<M-1>(in,grids);
    auto t1 = std::get<N-1-M>(grids)[in[N-1-M]];
    std::get<N-1-M>(out) = t1;
    return out;
}

template <typename ...GridTypes>
inline typename grid_tuple_traits<std::tuple<GridTypes...>>::indices 
    grid_tuple_traits<std::tuple<GridTypes...>>::get_indices(point_tuple in, const grid_tuple_type& grids)
{
    return get_indices<N-1>(in,grids);
}

template <typename ...GridTypes>
template <int M, typename std::enable_if<M ==0, bool>::type>
inline typename grid_tuple_traits<std::tuple<GridTypes...>>::indices 
    grid_tuple_traits<std::tuple<GridTypes...>>::get_indices(point_tuple in, const grid_tuple_type& grids)
{
    indices out;
    auto t1 = std::get<N-1>(in);
    #ifndef NDEBUG
    {
        if (std::get<N-1>(grids).size()<=t1.index_) throw std::logic_error("index mismatch");
        if (std::get<N-1>(grids)[t1.index_].index_ != t1.index_) { throw std::logic_error("index mismatch"); };
    }
    #endif
    out[N-1]=t1.index_;
    return out;
}

template <typename ...GridTypes>
template <int M, typename std::enable_if<M >= 1, bool>::type >
inline typename grid_tuple_traits<std::tuple<GridTypes...>>::indices 
    grid_tuple_traits<std::tuple<GridTypes...>>::get_indices(point_tuple in, const grid_tuple_type& grids)
{
    auto out = get_indices<M-1>(in,grids);
    auto t1 = std::get<N-1-M>(in);
    #ifndef NDEBUG
    if (std::get<N-1-M>(grids).size()<=t1.index_) throw std::logic_error("index mismatch");
    if (std::get<N-1-M>(grids)[t1.index_].index_ != t1.index_) { throw std::logic_error("index mismatch"); };
    #endif
    out[N-1-M]=t1.index_;
    return out;
}

//impl

template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...GridTypes, typename ...ArgTypes>
struct GridArgTypeExtractor<ValueType, T<GridType1, GridTypes...>, ArgTypes...> : 
GridArgTypeExtractor<ValueType, T<GridTypes...>, ArgTypes..., typename GridType1::point::value_type>
{
};

template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...ArgTypes>
struct GridArgTypeExtractor<ValueType, T<GridType1>, ArgTypes...> {
    typedef std::function<ValueType(ArgTypes...,typename GridType1::point::value_type)> type; 
    typedef std::function<ValueType(ArgTypes...,typename GridType1::point::value_type)> arg_type; 
    typedef std::tuple<ArgTypes...,typename GridType1::point::value_type> arg_tuple;
};


template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...GridTypes, typename ...ArgTypes>
struct GridPointExtractor<ValueType, T<GridType1, GridTypes...>, ArgTypes...> : 
GridPointExtractor<ValueType, T<GridTypes...>, ArgTypes...,typename GridType1::point>
{
};

template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...ArgTypes>
struct GridPointExtractor<ValueType, T<GridType1>, ArgTypes...> {
 //   typedef std::function<ValueType(ArgTypes...,typename GridType1::point)> point_type; 
    typedef std::tuple<ArgTypes...,typename GridType1::point> arg_tuple;
};


/* A tool to generate an array of grid sizes from a given tuple of grids. */
template <size_t N> 
struct GetGridSizes {
    template <typename ... GridType>
    static inline std::array<size_t,sizeof...(GridType)> TupleSizeToArray( const std::tuple<GridType...>& in ) {
        static_assert(N>1,"!");
        auto out = GetGridSizes<N-1>::template TupleSizeToArray<GridType...>( in );
        std::get<N-1>(out) = std::get<N-1>(in).size();
        return out;
    };
};

template <> 
struct GetGridSizes<1> {
    template <typename... GridType>
    static inline std::array<size_t, sizeof...(GridType)> TupleSizeToArray( const std::tuple<GridType...>& in ) {
        std::array<size_t,sizeof...(GridType)> out;
        std::get<0>(out) = std::get<0>(in).size();
        return out;
    }
};

/* A tool to convert a tuple of grid points to an array of indices */ 
template <typename ... Args> struct PointTupleToIndex {

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

} // end of namespace extra
} // end of namespace gftools

