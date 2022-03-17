#pragma once

#include <utility>
#include <array>
#include <cstddef>

//#include "tuple_tools.hpp" // only for debug

namespace gftools { 

/** eval_expression is an unary tree, representing operator[] of a multidimensional object V, 
  * which has a []...[] method of depth N. The evaluation takes place in 3 walks (forward-back-forward) in the tree,
  * using operator[], apply and get methods. operator[] returns a child eval_expression of depth D=N-1. 
  * When the depth D reaches 1, child calls apply(indices) method of a parent, passing the index of [] operator to it. 
  * Going back to the the root of the tree (N=D) apply method collects an array of indices. 
  * After that this array is passed to the static method get(data,indices) which calls operator [] of 
  * V recursively until D=1, gaining in the end the value of V at requested indices.
*/
template <typename V, typename VType, size_t D, size_t N = D> struct eval_expression 
{
    /// Typedef for underlying expression.
    typedef eval_expression<eval_expression<V,VType,D,N>, VType, D-1,N> under_type;
    /// Typedef of underylung value type of V
    typedef VType value_type;

    /// Construct from data object or eval_expression parent
    eval_expression(V v):v_(std::forward<V>(v)){};
    /// Returns a child - underlying expression of the order of D-1
    under_type operator[](int x) { i_ = x; return under_type(*this); }
    /// Combine indices from [] operator using apply method
    template <size_t M = N>
        typename std::enable_if<(M!=D), VType&>::type apply(std::array<int, N>&& indices){
            indices[N-D] = i_; return v_.apply(std::move(indices)); }
    /// Call [] operator of V using combined indices from eval_expression::apply method
    template <size_t M = N>
        typename std::enable_if<(M==D), VType&>::type apply(std::array<int, N>&& indices){
            indices[0] = i_; return this->get(v_, indices);
        };

    /// Pass indices, collected during first walk through the tree with [] and apply methods to V and use 
    template <typename V2>
    static VType& get(V2&& d, std::array<int, N> indices) { return under_type::get(d[indices[N-D]], indices); }

    /// Data or parent expression
    V v_; 
    /// Index of operator[]
    int i_;
};

template <typename V, typename VType, size_t N> struct eval_expression<V,VType,1,N>
{
    /// Construct from data object or eval_expression parent
    eval_expression(V v):v_(std::forward<V>(v)){};
    /// return value of v[x] for 1d objects.
    template <size_t M = N>
        typename std::enable_if<(M==1), VType&>::type operator[] (int x){return v_[x];}
    /// pass the index of operator[] to parent apply method if depth != 1
    template <size_t M = N>
        typename std::enable_if<(M!=1), VType&>::type operator[] (int x){std::array<int, N> ind; ind[N-1] = x; return v_.apply(std::move(ind));}

    /// Use collected indices to extract value from V
    template <typename V2>
        static VType& get(V2&& d, std::array<int, N> indices) { return d[indices[N-1]]; }

    /// Data or parent expression
    V v_;
};

};
