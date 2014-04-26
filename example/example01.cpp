
/* Here's a small example of the code, written with usage of some basic
 * objects of the library. It covers initialization of the grids and objects,
 * defined on the grid. Here's what's done:
 */

// This includes the grid of real values.
#include "real_grid.hpp"
// This includes the definition of a container class - a tool to store data.
#include "container.hpp"
// This include the grid_object class - an object to store data, defined on a grid.
#include "grid_object.hpp"

/* Some more basic definitions, that are used:
 * double is double
 * ComplexType is complex<double>
 * INFO( something ) is the same as std::cout << something << std::endl;
 * DEBUG( something ) is the same as INFO, but it also prints the number of the
 * line in the code, where it is called from.
 * ERROR( something ) is the same as INFO but is used for outputting the errors,
 * it sends the text to the error stream (stderr).
 */

// This is a shorthand for not writing GFTools::real_grid, GFTools::grid_object etc
using namespace gftools;

int main()
{
    // Let's start with grids.
    // This generates a uniform grid of 20 real numbers from -5 to 5
    real_grid grid1(-5,5,20);
    // And prints it
    INFO(grid1); // same as std::cout << grid1 << std::endl
    
    // New c++ feature - store a function
    std::function<double(int)> f1;
    /* The right hand side of next line is called a "lambda" construction.
     * We construct some function of integer and assign it to f1. */
    f1 = [](int n){return std::pow(3,n);};
    // f1 is now 3^n. Let's now generate a non-uniform grid and print it.
    real_grid grid2(0,10,f1);
    INFO(grid2)
    
    // Another option is to generate a set of values and generate a grid.
    std::vector<double> v1;
    for (double a=-3;a<=3;a+=0.01) v1.push_back(a);
    real_grid grid3(v1);
    INFO(grid3);
    
    // Now let's make an integration over a grid
    // Let's create another function
    std::function<double(double)> sinF;
    sinF = [](double x){return sin(x);};
    // sinF is now sin(x). We integrate an even function and generate zero.
    INFO(grid1.integrate(sinF));
    
    // Ok, that's all for grids.
    // Now, let's define objects on a grid.
    // First, make a shorthand
    typedef grid_object<double,real_grid> GFRetarded;
    // Now in order to construct an object on a grid, we need to provide a grid.
    // Let's use the one, already created.
    GFRetarded G1(grid1);
    // After construction, grid_object stores only zeros.
    // Let's fill it with something more interesting, say sin(x).
    G1.fill(sinF);
    INFO(G1)
    // And now let's integrate it. It should give zero once again.
    INFO(grid1.integrate(G1));
    
    // In order to have an access to the data points, one can just use [] operator.
    // Notice the values are counted from zero.
    INFO(G1[2]);
    
    //grid_objects also support mathematical operations. Let's create another
    //grid_object first.
    GFRetarded G2(grid1);
    std::function<double(double)> cosF;
    cosF = [](double x){return cos(x);};
    G2.fill(cosF);
    INFO(G2);
    // Now let's make some operations. sin^2(x) + cos^2(x) is trivial.
    GFRetarded G3 = G2*G2+G1*G1;
    INFO(G3);
    
    // Some advanced features. Grids allow to interpolate values.

#warning fix operator() of gridobjects
// fixme 
/*
    INFO(grid1.eval(G1,4.47));
    INFO(G1(4.47));
    INFO(grid1.eval(G1,0.0));
    INFO(G1(0.0));
    
    //grid_objects allow for the extrapolation of the values, if they are out 
    // of grid bounds. For this we need to assign the extrapolation
    // function. 
    G1.tail() = sinF;
    INFO(G1(2*PI));
    G1.savetxt("G1.dat");
*/

    grid_object<double,real_grid,real_grid> R1(std::make_tuple(grid1,grid1));
    std::function<double(double,double)> cos2;
    cos2 = [](double x,double y){return cos(x)*cos(y);};
    R1.fill(cos2);
    R1.savetxt("G2.dat");
}
