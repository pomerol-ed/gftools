#include <numeric>

#include "defaults.hpp"
#include "matsubara_grid.hpp"
#include "container.hpp"
#include "grid_object.hpp"
#include "enum_grid.hpp"

#include <iostream>
#include <ctime>
#include <array>

using namespace gftools;

int main()
{
    typedef grid_object<complex_type,fmatsubara_grid> GF;
    fmatsubara_grid n1(0,2,10);
    fmatsubara_grid n2(0,5,20);
    auto a1 = std::make_tuple(n1,n2);

    GF D1(n2);
    D1.getData()[3]=3.4;
    DEBUG(D1);
    DEBUG(D1(FMatsubara(3,20)));
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> D2(std::make_tuple(n1,n2));
    
    D2.getData()[0][1]=4.0;
    D2.getData()[1][2]=3.1;
    DEBUG(D2);
    D1.savetxt("D1.dat");
    D2.savetxt("D2.dat");

    auto& C1 = D2.getData();
    DEBUG(C1[0]);
    //decltype (C1[0]) x(std::make_tuple(3));
    container<complex_type,1> C2(C1[0]);
    DEBUG(C2[1]);
    DEBUG(n2.evaluate(C2, FMatsubara(1,20)));
    DEBUG(n1.evaluate(C1, FMatsubara(1,10)));
    auto C22 = n1.evaluate(C1, FMatsubara(0,10));
    DEBUG(C22+C2-C22*2.0);
    //DEBUG(n2.getValue(std::forward<Container<complex_type,1>>(n1.getValue(C1, FMatsubara(1,10))), FMatsubara(1,20)));

    DEBUG(D2(FMatsubara(1,10),FMatsubara(2,20)));

    std::function<complex_type(complex_type, complex_type)> f1 = [](const complex_type &a1, const complex_type &a2){return 1.0/(a1+a2);};
    DEBUG(f1(1.0,2.0));
    
    
    D2.fill(f1);
    DEBUG(D2);
    
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid, fmatsubara_grid> D3(std::make_tuple(n1,n2,n2));
    std::function<complex_type(complex_type, complex_type, complex_type)>  f2 = [](const complex_type &a1, const complex_type &a2, const complex_type &a3){return a3/(a1+a2);};
    std::function<complex_type(complex_type, complex_type, complex_type)>  f3 = [](const complex_type &a1, const complex_type &a2, const complex_type &a3){return a3*a2/(a1);};
    D3.fill(f2);
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid, fmatsubara_grid> D4(std::make_tuple(n1,n2,n2));
    D4.fill(f3);
    DEBUG(D4);
    DEBUG(D3);
    DEBUG(D3+D4);

    // Tails test.
    int n_freq = 6;
    real_type beta=10.0;
    fmatsubara_grid fgrid(-n_freq,n_freq,beta);
    real_type t=3.0;
    std::function<complex_type(complex_type)> f11, f21;
    f11 = [t](complex_type w) -> complex_type {return t*t/w;};
    f21 = [t](complex_type w) -> complex_type {return t;};
    grid_object<complex_type,fmatsubara_grid> D21(fgrid);
    grid_object<complex_type,fmatsubara_grid> D31(fgrid);
    D21 = t;
    DEBUG(D21);
    DEBUG(D21(FMatsubara(n_freq+1,beta)));

    D21.fill(f11);
    DEBUG(D21);
    DEBUG(D21(FMatsubara(n_freq+1,beta)));

    
    D21 = 2.0;
    DEBUG(D21(FMatsubara(n_freq+1,beta)));
    D21+=-1.0;
    DEBUG(D21(FMatsubara(n_freq+1,beta)));
    D31 = 3.0;
    DEBUG(D31(FMatsubara(n_freq+1,beta)));

    auto D41 = D21+D31;
    DEBUG(D41);
    DEBUG(D41(FMatsubara(n_freq+1,beta)));

    grid_object<real_type, enum_grid> g1(enum_grid(0,10));
    grid_object<real_type, enum_grid> g2(enum_grid(0,10));
    
    g1.fill(grid_object<real_type, enum_grid>::FunctionType([](int x){return x/2.0;}));
    DEBUG(g1);
    DEBUG(g1(3));
    DEBUG(g1(15));

    g1.savetxt("g1.dat");
    g2.loadtxt("g1.dat");
    return EXIT_SUCCESS;
}
