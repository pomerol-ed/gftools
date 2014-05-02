gftools is a set of tools to work with numeric physics-condmat problems. It provides three basic classes of objects:

- grids. 
There exists a set of implementation of grids for common purposes:
    - `matsubara_grid<b>`. A grid of Matsubara frequencies. "b" corresponds to a choice of fermionic(1) or bosonic(0) grids.
    - `real_grid `. A grid of floats, e.g. frequencies. 
    - `kmesh`. An equidistant periodic set of points (for sampling Brillouine zones).
    - `enum_grid`. A grid of integers.
- `container<N, ValueType>`. A recursive multidimensional container of depth <N> to store ValueType objects.
- `grid_object<Grids..., ValueType>`. Green's function. In a broader scope - a multidimensional object that is defined on a set of grids (Grids...) that contains ValueType values.
Examples :
    - `grid_object<ComplexType,matsubara_grid<1>>` - Green's function in imaginary frequencies
    - `grid_object<ComplexType,matsubara_grid<1>,kmesh>` - Green's function in imaginary frequencies and reciprocal 1d space.
    - `grid_object<ComplexType,real_grid>` - Retarded or Advanced Green's function, etc...

Look for examples in "example/" directory for detailed information.

The code is provided as a header-only library with a set of examples. Requires C++11, Boost headers and Eigen >=3.1 library.
