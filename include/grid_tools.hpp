#pragma once

#include "grid_base.hpp"
#include "tuple_tools.hpp"

namespace gftools {

namespace tools { 

/** A tool to extract ArgTypes(types of numbers) from a tuple of grids. */
template <typename ValueType, typename ... > struct GridArgTypeExtractor;
/** A tool to extract point types from a tuple of grids. */
template <typename ValueType, typename ... > struct GridPointExtractor;

/** A tuple of grids type traits. Allows to deduce underlying types of grids (points, values) 
    and convert between tuples of points, values and arrays of indices. */
template <typename> struct grid_tuple_traits;

template <typename ... GridTypes>
struct grid_tuple_traits<std::tuple<GridTypes...> >
{
    constexpr static size_t N = sizeof...(GridTypes);
    typedef std::tuple<GridTypes...> grid_tuple_type;
    /// Static sequence N-1, N-2, ... 1,0
    typedef typename tuple_tools::extra::index_gen<N>::type index_gen; 

    /// A typedef for a tuple of grid points.
    typedef typename GridPointExtractor<std::true_type, std::tuple<GridTypes...> >::arg_tuple point_tuple;
    /// A typedef for a tuple of grid point values.
    typedef typename GridArgTypeExtractor<std::true_type, std::tuple<GridTypes...> >::arg_tuple arg_tuple;
    /// A typedef for a set of indices. 
    typedef std::array<size_t, N> indices;

    /// Convert indices to points.
    static point_tuple get_points(indices in, const grid_tuple_type& grids);
    /// Convert indices to args.
    static arg_tuple get_args(indices in, const grid_tuple_type& grids);
    /// Convert points to indices.
    static indices get_indices(point_tuple in, const grid_tuple_type& grids);
    /// Get an array of dimensions of grids.
    static indices get_dimensions(const grid_tuple_type& grids) {return dims_(index_gen(), grids); }
    /// Get a product of all dimensions of gris.
    static int get_total_size(const grid_tuple_type& grids){auto d = get_dimensions(grids); return std::accumulate(d.begin(),d.end(),1,std::multiplies<int>()); }
    //template <class F> evaluate(F &&f, arg_tuple in, const grid_tuple_type& grids); 

protected:
    template <int...S> 
        static point_tuple get_points_(tuple_tools::extra::arg_seq<S...>, indices in, const grid_tuple_type& grids);
    template <int...S> 
        static arg_tuple get_args_(tuple_tools::extra::arg_seq<S...>, indices in, const grid_tuple_type& grids);
    template <int...S> 
        static indices get_indices_(tuple_tools::extra::arg_seq<S...>, point_tuple in, const grid_tuple_type& grids);
    template <int...S> 
        static indices dims_(tuple_tools::extra::arg_seq<S...>, const grid_tuple_type& grids) {return {{ (std::get<S>(grids).size())... }};};
};

template <typename ...GridTypes>
template <int...S> 
inline typename grid_tuple_traits<std::tuple<GridTypes...>>::arg_tuple 
    grid_tuple_traits<std::tuple<GridTypes...>>::get_args_(typename tuple_tools::extra::arg_seq<S...>, indices in, const grid_tuple_type& grids)
{
    #ifndef NDEBUG
    if ( indices({{ (std::get<S>(grids).size())... }}) <= indices({{ (std::get<S>(in))... }})) 
        throw std::logic_error("indices are out of grid bounds");
    #endif
    return std::make_tuple((std::get<S>(grids)[std::get<S>(in)].val_)...);
}

template <typename ...GridTypes>
inline typename grid_tuple_traits<std::tuple<GridTypes...>>::arg_tuple 
    grid_tuple_traits<std::tuple<GridTypes...>>::get_args(indices in, const grid_tuple_type& grids)
{
    return get_args_(index_gen(),in,grids);
}

template <typename ...GridTypes>
template <int...S> 
inline typename grid_tuple_traits<std::tuple<GridTypes...>>::point_tuple 
    grid_tuple_traits<std::tuple<GridTypes...>>::get_points_(typename tuple_tools::extra::arg_seq<S...>, indices in, const grid_tuple_type& grids)
{
    #ifndef NDEBUG
    if ( indices({{ (std::get<S>(grids).size())... }}) <= indices({{ (std::get<S>(in))... }})) 
        throw std::logic_error("indices are out of grid bounds");
    #endif
    return std::make_tuple((std::get<S>(grids)[std::get<S>(in)])...);
}

template <typename ...GridTypes>
inline typename grid_tuple_traits<std::tuple<GridTypes...>>::point_tuple 
    grid_tuple_traits<std::tuple<GridTypes...>>::get_points(indices in, const grid_tuple_type& grids)
{
    return get_points_(index_gen(),in,grids);
}

template <typename ...GridTypes>
template <int...S> 
inline typename grid_tuple_traits<std::tuple<GridTypes...>>::indices 
    grid_tuple_traits<std::tuple<GridTypes...>>::get_indices_(typename tuple_tools::extra::arg_seq<S...>, point_tuple in, const grid_tuple_type& grids)
{
    #ifndef NDEBUG
    if ( indices({{ (std::get<S>(grids).size())... }}) <= indices({{ (std::get<S>(in).index_)... }})) 
        throw std::logic_error("indices are out of grid bounds");
    if ( indices({{ (std::get<S>(grids)[std::get<S>(in).index_].index_) ... }}) != indices({{ (std::get<S>(in).index_)... }})) 
        throw std::logic_error("index mismatch");
    #endif
    return {{ (std::get<S>(in).index_)... }};
};

template <typename ...GridTypes>
inline typename grid_tuple_traits<std::tuple<GridTypes...>>::indices 
    grid_tuple_traits<std::tuple<GridTypes...>>::get_indices(point_tuple in, const grid_tuple_type& grids)
{
    return get_indices_(index_gen(),in,grids);
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

// obsolete

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

