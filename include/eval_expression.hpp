#pragma once

#include "grid_tools.hpp"

namespace gftools { 

// eval some object on a grid

template <typename Node, typename, typename Evaluable = Node> struct grid_eval_expr; 

template <typename Node, typename ...GridTypes, typename Evaluable> 
struct grid_eval_expr<Node, std::tuple<GridTypes...>, Evaluable>
{ 
    constexpr static size_t D = sizeof...(GridTypes);
    typedef tools::grid_tuple_traits<std::tuple<typename std::remove_reference<GridTypes>::type...>> trs;
    typedef typename std::tuple<typename std::add_const<typename std::add_lvalue_reference<GridTypes>::type>::type...> grid_constref_tuple; 
    typedef typename trs::grid_tuple_type grid_tuple_type;
    typedef typename trs::arg_tuple arg_tuple;
    typedef typename std::remove_reference<decltype(std::declval<Evaluable>()[0])>::type data_type;
    typedef grid_eval_expr< 
        grid_eval_expr<Node, std::tuple<GridTypes...>> const&, 
        decltype(tuple_tools::tuple_tail(std::declval<grid_constref_tuple>())),
        data_type
    > child_t;

    Node f_;
    grid_constref_tuple grids_;

    grid_eval_expr() = delete;
    grid_eval_expr(Node F, grid_constref_tuple grids):
        f_(std::forward<Node>(F)),
        grids_(grids)
    {
    DEBUG(D);
    }

    grid_eval_expr(grid_eval_expr const&) = delete;
    grid_eval_expr& operator= (grid_eval_expr const&) = delete;

    grid_eval_expr(grid_eval_expr&&) = default;
    grid_eval_expr& operator= (grid_eval_expr&&) = default;

    //template <bool Root = (std::is_same<Node,Evaluable>::value), typename std::enable_if<!Root>::type>

    template <typename T = Node, typename std::enable_if<!std::is_same<T,Evaluable>::value>::type>
    auto apply (std::array<size_t, D-1> p, size_t i) -> decltype(f_.apply({{i}}, 0)) {}

    template <typename T = Node, typename std::enable_if<std::is_same<Node,Evaluable>::value>::type>
    auto apply (std::array<size_t, D-1> p, size_t i) const -> int {} ;//decltype(f_.apply({{i}}, 0)) {}

    template <int D1 = D>
    typename std::enable_if<(D1>1),child_t>::type
    operator[] (size_t index) const { 
        return child_t(*this,tuple_tools::tuple_tail(grids_)); 
    }

    template <int D1 = D>
    typename std::enable_if<(D1==1), data_type>::type
    operator[] (size_t index) const { 
        return apply({{}}, index);
    }
    

    //auto eval(arg_tuple_type in) -> child_t const;

};

};
