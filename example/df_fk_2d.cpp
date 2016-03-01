/** 
 * This is a simple unoptimized code to evaluate some low-order dual diagrams in FK model
 * input : 
 *  --order : order of the diagram
 *  --vertex_file : full impurity vertex. default : "gamma4.dat"
 *  --gd0_file : bare dual Green's function. default : "gd0_k.dat"
 * output:
 *  - dual self-energy "sigma_wk.dat"
 *  - cut of dual self-energy at first Matsubara freq
 *  - k-dependence of dual bubbles (summer over Matsubara freqs). 
 */

#include <boost/program_options.hpp>
#include <gftools.hpp>
#include <Eigen/Core>

namespace po = boost::program_options;
using namespace gftools;

int main(int argc, char *argv[])
{
    // parse command line options
    po::options_description desc("FK diagrams evaluator"); 
    desc.add_options()
        ("order,n", po::value<int>()->default_value(2), "order of diagrams")
        ("beta", po::value<double>()->default_value(1), "inverse temperature")
        ("U", po::value<double>()->default_value(1), "inverse temperature")
        ("help", "produce help message");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) { std::cout << desc << std::endl; exit(0); }
    int diagram_order = vm["order"].as<int>(); 
    // we will use diagram_order as amount of iterations in vertices. 1 iteration = 2nd order, etc. 
    diagram_order-=(diagram_order>0);
    double beta = vm["beta"].as<double>(); 
    double U = vm["U"].as<double>(); 

    std::cout << "!" << std::endl;

    static constexpr int NDim = 2; // work in 2 dimensions
    int kpts = 16;

    typedef grid_object<std::complex<double>, fmatsubara_grid, fmatsubara_grid> vertex_type;
    typedef grid_object<std::complex<double>, fmatsubara_grid> gw_type;
    typedef grid_object<std::complex<double>, fmatsubara_grid, kmesh, kmesh> gk_type;
    typedef grid_object<std::complex<double>, kmesh, kmesh> disp_type;
    typedef Eigen::MatrixXcd matrix_type; 

    double T = 1.0/beta;
    std::cout << "beta = " << beta << std::endl;
    std::cout << "T = " << T << std::endl;

    // Prepare evaluation grids
    fmatsubara_grid fgrid(-20,20,beta);

    // Hybridization function (atomic limit)
    gw_type delta(fgrid);
    delta.fill([&](std::complex<double> w){return double(2*NDim) / w; });
    // Local gf (atomic limit)
    gw_type gw(fgrid);
    gw.fill([&](std::complex<double> w){return 0.5/(w - U/2.0) + 0.5/(w + U/2.0);});

    // Local vertex (atomic limit)
    vertex_type gamma4(std::make_tuple(fgrid, fgrid));
    typename gw_type::function_type Lambda = [U](std::complex<double> w){return 1. - U*U/4./w/w;};
    gamma4.fill([&](std::complex<double> w1, std::complex<double> w2){ 
            return beta * U * U / 4.0 * double(1 - tools::is_float_equal(w1, w2) ) * (1. - U*U/4./w1/w1) * (1. - U*U/4./w2/w2); 
        });

    // k-dependent input
    // define a mesh in Brilloin zone (BZ)
    kmesh kgrid(kpts);
    int totalkpts = std::pow(kpts, NDim);
    std::cout << "kmesh : " << kpts << " points " << std::endl; 
    std::cout << "total pts in BZ = " << totalkpts << std::endl; 

    // lattice dispersion
    disp_type eps(std::forward_as_tuple(kgrid, kgrid));
    eps.fill([&](double kx, double ky){return -2*(cos(kx) + cos(ky));});

    // lattice gf in dmft
    gk_type glat_dmft(std::forward_as_tuple(fgrid, kgrid, kgrid));
    glat_dmft.fill([&](fmatsubara_grid::point w, kmesh::point kx, kmesh::point ky){return 1.0 / (1.0 / gw(w) + delta(w) - eps(kx, ky));});

    // bare dual df
    gk_type gd0(glat_dmft.grids());
    gd0.fill([&](fmatsubara_grid::point w, kmesh::point kx, kmesh::point ky){return glat_dmft(w, kx, ky) - gw(w);});
    
    // define full vertex (will be evaluated below)
    vertex_type full_vertex(gamma4);
    // make a matrix from the vertex
    matrix_type gamma4_matrix = gamma4.data().as_matrix();
    matrix_type full_vertex_matrix(gamma4_matrix);

    // check that we are consistent - 2 frequency grids are the same
    if (gd0.template grid<0>() != fgrid) throw std::logic_error("matsubara grid mismatch");
    
    std::cout << "kmesh : " << kgrid << std::endl;
    // Initalize dual self-energy on the same grids as gd0
    gk_type sigma_dual(gd0.grids());
    sigma_dual = 0.0;
    disp_type bare_bubbles(kgrid,kgrid);

    // now loop through the BZ (no irreducible part optimization)
    for (kmesh::point q1 : kgrid.points()) { 
        for (kmesh::point q2 : kgrid.points()) { 
            std::cout << "[" << q1.index()*kpts + q2.index()+1 << "/" << totalkpts <<  "]; q = {" <<  q1.value() << " " << q2.value() << "}" << std::endl;
            // Evaluate -T \sum_k G_{w,k} G_{w,k+q}
            // Shift dual g in k-space and make no frequency shift. 
            // note : this operation is typically optimized via an fft. 
            gk_type gd0_shift = gd0.shift(std::make_tuple(0.0,q1,q2)); 
            // obtain a bubble 
            gk_type bubble_wk = -T * gd0 * gd0_shift;
            // perform sum over k    
            gw_type dual_bubble(fgrid);
            for (auto w : fgrid.points()) { dual_bubble[w] = bubble_wk[w].sum() / double(totalkpts); }
            // save the bubble for output 
            bare_bubbles(q1, q2) = dual_bubble.sum();
            // construct a diagonal matrix (in frequency space from the bubble)
            matrix_type dual_bubble_matrix = dual_bubble.data().as_diagonal_matrix(); 
            // refresh vertex
            full_vertex_matrix = gamma4_matrix;
            // IMPORTANT : Evaluate dual diagram
            for (int n=0; n<diagram_order; n++) 
                full_vertex_matrix= gamma4_matrix * dual_bubble_matrix * full_vertex_matrix; 
            // update self-energy
            for (auto w : fgrid.points())
                sigma_dual[w.index()] += T* full_vertex_matrix(w.index(), w.index()) * gd0_shift[w.index()] / double(totalkpts);
        }
    }
    // output
    // save sigma
    sigma_dual.savetxt("sigma_wk.dat");
    // save sigma at first matsubara
    auto w0 = fgrid.find_nearest(I*PI/beta);
    disp_type sigma_w0(std::forward_as_tuple(kgrid,kgrid),sigma_dual[w0]);
    sigma_w0.savetxt("sigma_w0.dat");
    // save bare bubbles
    bare_bubbles.savetxt("db0.dat");
}

