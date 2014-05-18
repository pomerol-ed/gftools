#pragma once

#include <utility>

namespace gftools { 

template <typename L, typename Op, typename R>
struct math_expr { 

    math_expr(L l, R r):l_(std::forward<L>(l)),r_(std::forward<R>(r)){}
    math_expr(math_expr const&) = delete;
    math_expr& operator=(math_expr const&) = delete;
    math_expr(math_expr&& rhs):l_(std::forward<L>(rhs.l_)),r_(std::forward<R>(rhs.r_)){}
    math_expr& operator=(math_expr &&) = default;

    template <typename RHS> 
    using result_type = math_expr<math_expr<L,Op,R>,Op,RHS&&>; 

    L l_;
    R r_;
};

} // end of namespace gftools
