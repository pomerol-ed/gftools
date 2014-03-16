#include <numeric>

#include "grid_base.hpp"

using namespace gftools;

struct data1 { double x, y; data1(double z):x(z),y(z*2){}; operator double() const {return x;} };
bool operator< (data1 x, data1 y){return x.x < y.x;};
std::ostream& operator<<(std::ostream& out, data1 d){out << "(" << d.x << " " << d.y << ")"; return out;};

struct my_grid : public grid_base<data1, my_grid> {
    typedef grid_base<data1, my_grid> base;
    using base::vals_;
    my_grid() = default;
    my_grid(std::vector<value_type> v):base(v){};
    
    std::tuple <bool, size_t, real_type> find (data1 in) const {
        auto out = std::lower_bound (vals_.begin(), vals_.end(), in);
        size_t i = size_t(out-vals_.begin());
        i--;
        if (i==vals_.size()-1) return std::make_tuple(1,i,1.0);
        return std::make_tuple (1,i,1.0);
    }
};

int main()
{

    my_grid g1;
    
    my_grid g2 ( { data1(1), data1(2), data1(5) });
    std::cout << g2 << std::endl;
    for (auto p : g2.get_points()) std::cout << p << " "; std::cout << std::endl;

    typedef my_grid::point point;
    point p1 = g2[2];
    std::cout << p1 << std::endl;

    point p2 = g2.find_closest(data1(3));
    std::cout << p2 << " == " << g2[1] << " = " << std::boolalpha << (p2 == g2[1]) <<  std::endl;
    if (!(p2 == g2[1])) return EXIT_FAILURE;
/*
    FMatsubaraGrid n1(-100,100,10);
    FMatsubaraGrid n2(0,32,20);

    typename FMatsubaraGrid::point ppp = n1[0];
    std::cout << ppp << std::endl;    


    KMesh k1(32);
    std::function<RealType(RealType)> F1=[](RealType x) {return cos(x);};
    RealType outF1 = k1.integrate (F1);
    INFO(outF1);
    if (std::abs(outF1)>1e-8) return EXIT_FAILURE;

    std::function<ComplexType(ComplexType)> F2=[](ComplexType x) {return x;};
    ComplexType out = n1.integrate(F2);
    if (std::abs(out)>1e-8) return EXIT_FAILURE;
    INFO(out);
    DEBUG(n2);

    auto a1 = std::make_tuple(n1,n2);
    
    std::vector<int> x(n2.getSize());
    for (size_t i=0; i<x.size(); ++i) x[i]=i*i;

    INFO(std::get<1>(n2.find(FMatsubara(10,20))));
    if (std::get<1>(n2.find(FMatsubara(10,20)))!=10) return EXIT_FAILURE;
    INFO(std::get<1>(n2.find(FMatsubara(20,20))));
    if (std::get<1>(n2.find(FMatsubara(20,20)))!=20) return EXIT_FAILURE;
    INFO(std::get<1>(n2.find(FMatsubara(32,20))));
    if (std::get<1>(n2.find(FMatsubara(32,20)))!=32) return EXIT_FAILURE;
    INFO(std::get<1>(n2.find(FMatsubara(35,20))));

    INFO(n2.getValue(x,FMatsubara(10,20)));
    if (n2.getValue(x,FMatsubara(10,20))!=100) return EXIT_FAILURE;

    RealGrid grid1(-5,5,20);
    // And prints it
    INFO(grid1); // same as std::cout << grid1 << std::endl
    
    // New c++ feature - store a function
    std::function<RealType(int)> f1;
    f1 = [](int n){return std::pow(3,n);};
    // f1 is now 3^n. Let's now generate a non-uniform grid and print it.
    RealGrid grid2(0,10,f1);
    INFO(grid2)
    
    // Another option is to generate a set of values and generate a grid.
    std::vector<RealType> v1;
    for (RealType a=-3;a<=3;a+=0.01) v1.push_back(a);
    RealGrid grid3(v1);
    INFO(grid3);
    
    // Now let's make an integration over a grid
    // Let's create another function
    std::function<RealType(RealType)> sinF;
    sinF = [](RealType x){return sin(x);};
    // sinF is now sin(x). We integrate an even function and generate zero.
    INFO(grid1.integrate(sinF));
    
 

    EnumerateGrid g1(3,8);
    INFO(g1);
    EnumerateGrid::point p1 = g1[0];
    std::cout << p1 << std::endl;
*/
    return EXIT_SUCCESS;
}
