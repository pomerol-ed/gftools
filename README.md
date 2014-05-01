This is a set of tools to work with numeric condmat problems. It consists of three basic objects:

- grids. 
There exists a set of implementation of grids for common purposes:
    - matsubara_grid< b >. A grid of Matsubara frequencies. "b" corresponds to a choice of fermionic(1) or bosonic(0) grids.
    - real_grid . A grid of floats, e.g. frequencies. 
    - kmesh. An equidistant set of points from 0 to 2pi.
    - enum_grid. A grid of integers.
- container<N, ValueType>. A recursive multidimensional container of depth <N> to store ValueType objects.
- grid_object<Grids..., ValueType>. An multidimensional object that is defined on a set of grids, Grids... that contains ValueType values.
Example :
    - grid_object<ComplexType,matsubara_grid<1>> - Green's function in imaginary frequencies
    - grid_object<ComplexType,matsubara_grid<1>,KMesh> - Green's function in imaginary frequencies and reciprocal 1d space.
    - grid_object<ComplexType,RealGrid> - Retarded or Advanced Green's function, etc...

Look for examples in "example/" directory for detailed information.

The code is provided as a header-only library with a set of examples. Requires C++11 and >=Eigen3.1 library.
