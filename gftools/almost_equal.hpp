#pragma once

namespace gftools {

//almost equal comparison
///TODO: note duplication of functionality with the more fancy is_float_equal function in tools.
template <typename T1,typename T2>
    bool almost_equal(T1 t1, T2 t2, real_type tol = 10.*std::numeric_limits<real_type>::epsilon()){ 
  return std::abs(t1-t2)<tol; 
}

} // end namespace gftools

