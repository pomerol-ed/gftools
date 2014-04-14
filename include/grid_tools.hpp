#pragma once

#include "grid_base.hpp"

namespace gftools {

template <typename ValueType, typename ... > struct GridArgTypeExtractor;
template <typename ValueType, typename ... > struct GridPointExtractor;

/** A tool to generate a function of argtypes of grids. */
template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...GridTypes, typename ...ArgTypes>
struct GridArgTypeExtractor<ValueType, T<GridType1, GridTypes...>, ArgTypes...> : 
GridArgTypeExtractor<ValueType, T<GridTypes...>, ArgTypes..., typename GridType1::point::value_type>
{
};

template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...ArgTypes>
struct GridArgTypeExtractor<ValueType, T<GridType1>, ArgTypes...> {
    typedef std::function<ValueType(ArgTypes...,typename GridType1::point::value_type)> type; 
    typedef std::function<ValueType(ArgTypes...,typename GridType1::point::value_type)> arg_type; 
    typedef std::tuple<ArgTypes...,typename GridType1::point::value_type> arg_tuple_type;
    typedef tuple_tools::extra::tuple_caller<ValueType, std::tuple<ArgTypes..., typename GridType1::point::value_type>> arg_function_wrapper;
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
    typedef tuple_tools::extra::tuple_caller<ValueType, std::tuple<ArgTypes..., typename GridType1::point>> point_function_wrapper;
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

} // end of namespace gftools

