#pragma once

#include <tuple>
#include <array>
#include <vector>

#include "defaults.hpp"
#include "num_io.hpp"

//#include<memory>
//#include<utility>
//#include<functional>

namespace gftools {

namespace tuple_tools {

//
// Split a tuple - a tool from http://stackoverflow.com/questions/10626856/how-to-split-a-tuple. */
//

namespace extra {
    template <typename Target, typename Tuple, int N, bool end >
    struct split_tuple_struct
    {
        template < typename ... Args >
        static Target create(Tuple const& t, Args && ... args)
        {
            return split_tuple_struct<Target,Tuple, N+1, std::tuple_size<Tuple>::value == N+1>::create(t, std::forward<Args>(args)..., std::get<N>(t));
        }
    };

    template < typename Target, typename Tuple, int N >
    struct split_tuple_struct<Target,Tuple,N,true>
    {
        template < typename ... Args >
        static Target create(Tuple const& t, Args && ... args) { return Target(std::forward<Args>(args)...); }
    };
};

/// Get the tail of the tuple
template < typename Head, typename ... Tail >
inline std::tuple<Tail...> tuple_tail(std::tuple<Head,Tail...> const& tpl)
{
    return extra::split_tuple_struct<std::tuple<Tail...>, std::tuple<Head,Tail...>, 1, std::tuple_size<std::tuple<Head,Tail...>>::value == 1>::create(tpl);
}

//
// Print a tuple
//

namespace extra {
    template <class T> struct tuple_io;
    template <typename Arg1,typename ...ArgTypes> struct tuple_io<std::tuple<Arg1,ArgTypes...>> 
        { 
            typedef std::tuple<Arg1,ArgTypes...> tuple_type;
            typedef Arg1 head_type;
            typedef std::tuple<ArgTypes...> tail_type;

            static std::string print(tuple_type in){ 
                std::stringstream out;
                out << std::get<0>(in) << " " << std::flush; 
                tail_type a = tuple_tail(in); 
                std::string t1 = tuple_io<tail_type>::print(a); 
                out << t1;
                return out.str();
            };

            static std::string serialize(std::tuple<Arg1,ArgTypes...> in){ 
                std::string out; std::stringstream b; 
                auto v = std::get<0>(in);
                b << num_io<Arg1>(v);
                out = b.str();
                tail_type d = tuple_tail(in);
                std::string t_print = tuple_io<tail_type>::serialize(d); 
                out+=" "; out+=t_print;
                return out;
                };
            static std::tuple<Arg1,ArgTypes...> read (std::istream &in) {
                head_type out; 
                num_io<head_type> t1(out);
                in >> t1;
                return std::tuple_cat(std::make_tuple(t1()),tuple_io<tail_type>::read(in));
                }
        };
    template <typename ArgType> struct tuple_io<std::tuple<ArgType>> 
        { 
            static std::string print(std::tuple<ArgType> in ) {
                std::stringstream s1; s1 << std::get<0>(in); return s1.str();
                };
            static std::string serialize(std::tuple<ArgType> in){ 
                std::string out; std::stringstream b; 
                auto v = std::get<0>(in);
                b << num_io<decltype(v)>(v) << std::flush;
                return b.str();
            };
            static std::tuple<ArgType> read (std::istream &in) {
                ArgType out;
                num_io<ArgType> t1(out);
                in >> t1;
                return std::make_tuple(t1());
            };
        };
}; // end of namespace extra

/// Serialize a tuple using num_io
template <typename TT>
std::string serialize_tuple(TT t1){
    return std::forward<std::string>(extra::tuple_io<TT>::serialize(t1));};

/// Get a print of tuple via stream << on all members
template <typename TT>
std::string print_tuple(TT t1){
    return std::forward<std::string>(extra::tuple_io<TT>::print(t1));};

template <typename Arg, size_t D>
std::string print_array(const std::array<Arg,D>& t1){
    std::stringstream out; for (auto x:t1) out << x << " " << std::flush; return std::move<std::string>(out.str()); }

/// Read a tuple from a stream;
template<typename TT>
TT read_tuple(std::istream &in){return extra::tuple_io<TT>::read(in);};


// helpers for unfolding
namespace extra {
    template<int ...> struct arg_seq {};
    template<int N, int ...S> struct index_gen : index_gen<N-1, N-1, S...> {};
    template<int ...S> struct index_gen<0, S...>{ typedef arg_seq<S...> type; };
}; 

//
/// create a tuple or array from by repeating an argument
//

template <typename Arg, size_t D, typename... Extras> struct repeater : repeater<Arg,D-1,Arg,Extras...>  {  
    template <int ...S, typename ...ArgTypes> static std::array<Arg,D> get_arr_(extra::arg_seq<S...>, std::tuple<ArgTypes...> in) { return  {{ std::get<S>(in)... }}; }
    static std::array<Arg,D> get_array(const Arg& in){ return get_arr_(typename extra::index_gen<D>::type(), get_tuple(in)); };
    static typename repeater::tuple_type get_tuple(const Arg& in){ return std::tuple_cat(std::forward_as_tuple(in),repeater<Arg,D-1>::get_tuple(in)); }; 
};

template <typename Arg, typename ... Extras> struct repeater<Arg,1,Extras...> { 
    typedef std::tuple<Arg,Extras...> tuple_type;
    static std::array<Arg,1> get_array(const Arg& in){ return {{ in }}; };
    static tuple_type get_tuple(const Arg& in){ return std::forward_as_tuple(in); }
};



/** Caller for a function from the given tuple. */
//template <typename ReturnType, typename ...Args> struct __caller { 
//    std::tuple<Args...> _params;
//    std::function<ReturnType(Args...)> _f;
//    template<int ...S> ReturnType _callf(extra::arg_seq<S...>) { return _f(std::get<S>(_params)...); };
//    ReturnType call(){ return _callf(typename extra::index_gen<sizeof...(Args)>::type()); };
//};

namespace extra {
    template <typename, typename > struct tuple_caller;
    template <typename ReturnType, typename ...Args> struct tuple_caller<ReturnType,std::tuple<Args...>> { 
        typedef std::function<ReturnType(Args...)> function_type;
        typedef std::tuple<Args...> tuple_type;

        tuple_caller(const function_type& F, tuple_type in):f_(F),params_(in){};
        template<int ...S> ReturnType callf_(extra::arg_seq<S...>) { return f_(std::get<S>(params_)...); };
        ReturnType call(){ return callf_(typename extra::index_gen<sizeof...(Args)>::type()); };

        const function_type& f_;
        tuple_type params_;
    };
}


//
// unfold tuple into a function
//
template <typename Functor, typename ... Args>
typename std::result_of<Functor(Args...)>::type unfold_tuple(Functor F, std::tuple<Args...> t) {
    typedef typename std::result_of<Functor(Args...)>::type result_type; 
    return extra::tuple_caller<result_type,std::tuple<Args...>>(F, t).call();
};

namespace extra { 
template <typename Functor, typename V, int...S> 
    typename Functor::result_type unpack_vector_(Functor F, std::vector<V> v, extra::arg_seq<S...>) { return F(v[S]...); } 
template <typename V, int ...S> 
    typename repeater<typename V::value_type, sizeof...(S)>::tuple_type array_to_tuple_(V v, extra::arg_seq<S...>) { return std::make_tuple(v[S]...); }
}

///
/// Unfold vector into a function
///
template <typename Functor, typename Arg, int N>
decltype(unfold_tuple(std::declval<Functor>(), std::declval< typename repeater<Arg, N>::tuple_type >())) 
    unfold_vector(Functor F, std::vector<Arg> t) {
    return extra::unpack_vector_(F, t, typename extra::index_gen<N>::type()); 
};

/// Convert std::array<V,N> to tuple<V...>. Normally compilers should do it automatically, but GCC and intel fails, so here's a cast.
template <typename V> 
typename repeater<typename V::value_type, std::tuple_size<V>::value>::tuple_type array_to_tuple(V x) {
    return extra::array_to_tuple_(x, typename extra::index_gen<std::tuple_size<V>::value>::type()); 
}


} // end of namespace tuple_tools

} // end namespace gftools


