#include <numeric>

#include "grid_base.hpp"

using namespace gftools;

struct my_data { double x, y; my_data(double z):x(z),y(z*2){}; operator double() const {return x;} };
bool operator< (my_data x, my_data y){return x.x < y.x;};
std::ostream& operator<<(std::ostream& out, my_data d){out << "(" << d.x << " " << d.y << ")"; return out;};

struct my_grid : public grid_base<my_data, my_grid> {
    typedef grid_base<my_data, my_grid> base;
    using base::vals_;
    my_grid() = default;
    my_grid(std::vector<value_type> v):base(v){};
    
    std::tuple <bool, size_t, real_type> find (my_data in) const {
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
    
    my_grid g2 ( { my_data(1), my_data(2), my_data(5) });
    std::cout << g2 << std::endl;
    for (auto p : g2.points()) std::cout << p << " ";
    std::cout << std::endl;

    typedef my_grid::point point;
    point p1 = g2[2];
    std::cout << p1 << std::endl;

    point p2 = g2.find_nearest(my_data(3));
    std::cout << p2 << " == " << g2[1] << " = " << std::boolalpha << (p2 == g2[1]) <<  std::endl;
    if (!(p2 == g2[1])) return EXIT_FAILURE;

    std::cout << g2[1] << " == " << g2.points()[1] << std::endl;
    if (g2[1]!=g2.points()[1]) return EXIT_FAILURE;

    std::cout << "grid size = " << g2.size() << std::endl;
    if (g2.size() != 3) return EXIT_FAILURE;

    std::vector<int> data_y ( { 1, 4, 5} );
    bool success =  (g2.eval(data_y, g2[1]) == data_y[1]);
    std::cout << "y[" << g2[1] << "]=" << g2.eval(data_y, g2[1]) << std::endl;
    if (!success) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
