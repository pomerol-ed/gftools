#pragma once 

#include <gftools/container.hpp>
#include <fftw3.h>

namespace gftools {

template <size_t D, typename BC, typename std::enable_if<D==1, bool>::type=0> 
container<complex_type,D> run_fft (const container_base<complex_type,D,BC> &in, int direction)
{
    container<complex_type,D> out(in);
    fftw_plan p;
    std::array<size_t,1> shape = {{ out.shape()[0] }};
    p = fftw_plan_dft_1d(shape[0], 
                         reinterpret_cast<fftw_complex*>(out.data()),
                         reinterpret_cast<fftw_complex*>(out.data()),
                         direction, FFTW_ESTIMATE); 
    real_type norm=1.0*shape[0];
    fftw_execute(p);
    if (direction == FFTW_BACKWARD) out/=norm;
    return out;
}


template <size_t D, typename BC, typename std::enable_if<D==2, bool>::type=0> 
container<complex_type,D> run_fft (const container_base<complex_type,D,BC> &in, int direction)
{
    container<complex_type,D> out(in);
    fftw_plan p;
    std::array<size_t,2> shape = {{ out.shape()[0], out.shape()[1] }};
    p = fftw_plan_dft_2d(shape[0], shape[1],  
                         reinterpret_cast<fftw_complex*>(out.data()),
                         reinterpret_cast<fftw_complex*>(out.data()),
                         direction, FFTW_ESTIMATE); 
    real_type norm=1.0*shape[0]*shape[1];
    fftw_execute(p);
    if (direction == FFTW_BACKWARD) out/=norm;
    return out;
}

template <size_t D, typename BC, typename std::enable_if<D==3, bool>::type=0> 
container<complex_type,D> run_fft (const container_base<complex_type,D,BC> &in, int direction)
{
    container<complex_type,D> out(in);
    fftw_plan p;
    std::array<size_t,3> shape = {{ out.shape()[0], out.shape()[1], out.shape()[2] }};
    p = fftw_plan_dft_3d(shape[0], shape[1], shape[2], 
                         reinterpret_cast<fftw_complex*>(out.data()),
                         reinterpret_cast<fftw_complex*>(out.data()),
                         direction, FFTW_ESTIMATE); 
    real_type norm=1.0*shape[0]*shape[1]*shape[2];
    fftw_execute(p);
    if (direction == FFTW_BACKWARD) out/=norm;
    return out;
}

template <size_t D, typename BC, typename std::enable_if<D==4, bool>::type=0> 
container<complex_type,D> run_fft (const container_base<complex_type,D,BC> &in, int direction)
{
    container<complex_type,D> out(in);
    fftw_plan p;
    const std::array<int,4> shape = {{ static_cast<int>(out.shape()[0]), static_cast<int>(out.shape()[1]),
                                       static_cast<int>(out.shape()[2]), static_cast<int>(out.shape()[3])
                                    }};
    p = fftw_plan_dft(4, shape.data(),
                         reinterpret_cast<fftw_complex*>(out.data()),
                         reinterpret_cast<fftw_complex*>(out.data()),
                         direction, FFTW_ESTIMATE); 
    real_type norm=1.0*shape[0]*shape[1]*shape[2]*shape[3];
    fftw_execute(p);
    if (direction == FFTW_BACKWARD) out/=norm;
    return out;
}



template <size_t D, typename BC, typename std::enable_if<D>=5, bool>::type=0> 
container<complex_type,D> run_fft (const container_base<complex_type,D,BC> &in, int direction)
{
    ERROR("No FFT defined for D="<<D);
    return in; 
}

} // end of namespace gftools
