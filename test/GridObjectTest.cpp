#include <numeric>

#include "Defaults.hpp"
#include "MatsubaraGrid.hpp"
#include "Container.hpp"
#include "GridObject.hpp"

#include <iostream>
#include <ctime>
#include <array>

using namespace GFTools;

int main()
{
    Log.setDebugging(true);
    typedef GridObject<ComplexType,FMatsubaraGrid> GF;
    FMatsubaraGrid n1(0,2,10);
    FMatsubaraGrid n2(0,5,20);
    auto a1 = std::make_tuple(n1,n2);

    GF D1(n2);
    D1.getData()[3]=3.4;
    DEBUG(D1);
    DEBUG(D1(FMatsubara(3,20)));
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> D2(std::make_tuple(n1,n2));
    
    D2.getData()[0][1]=4.0;
    D2.getData()[1][2]=3.1;
    DEBUG(D2);

    auto& C1 = D2.getData();
    DEBUG(C1[0]);
    //decltype (C1[0]) x(std::make_tuple(3));
    Container<1,ComplexType> C2(C1[0]);
    DEBUG(C2[1]);
    DEBUG(n2.getValue(C2, FMatsubara(1,20)));
    DEBUG(n1.getValue(C1, FMatsubara(1,10)));
    auto C22 = n1.getValue(C1, FMatsubara(0,10));
    DEBUG(C22+C2-C22*2.0);
    DEBUG(n2.getValue(n1.getValue(C1, FMatsubara(1,10)), FMatsubara(1,20)));

    DEBUG(D2(FMatsubara(1,10),FMatsubara(2,20)));

    std::function<ComplexType(ComplexType, ComplexType)> f1 = [](const ComplexType &a1, const ComplexType &a2){return 1.0/(a1+a2);};
    DEBUG(f1(1.0,2.0));
    
    
    D2.fill(f1);
    DEBUG(D2);
    
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid, FMatsubaraGrid> D3(std::make_tuple(n1,n2,n2));
    std::function<ComplexType(ComplexType, ComplexType, ComplexType)>  f2 = [](const ComplexType &a1, const ComplexType &a2, const ComplexType &a3){return a3/(a1+a2);};
    std::function<ComplexType(ComplexType, ComplexType, ComplexType)>  f3 = [](const ComplexType &a1, const ComplexType &a2, const ComplexType &a3){return a3*a2/(a1);};
    D3.fill(f2);
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid, FMatsubaraGrid> D4(std::make_tuple(n1,n2,n2));
    D4.fill(f3);
    DEBUG(D4);
    DEBUG(D3);
    DEBUG(D3+D4);

    /** Tails test. */
    int n_freq = 6;
    RealType beta=10.0;
    FMatsubaraGrid fgrid(-n_freq,n_freq,beta);
    RealType t=3.0;
    std::function<ComplexType(ComplexType)> f11, f21;
    f11 = [t](ComplexType w) -> ComplexType {return t*t/w;};
    f21 = [t](ComplexType w) -> ComplexType {return t;};
    GridObject<ComplexType,FMatsubaraGrid> D21(fgrid);
    GridObject<ComplexType,FMatsubaraGrid> D31(fgrid);
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
 
    return EXIT_SUCCESS;
}
