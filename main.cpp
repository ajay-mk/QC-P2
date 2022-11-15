#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigenvalues>

// Libint Gaussian integrals library
#include <libint2.hpp>
#if !LIBINT2_CONSTEXPR_STATICS
#include <libint2/statics_definition.h>
#endif

// Typedefs
using real_t = libint2::scalar_type;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

struct params {
    std::string type;
    std::string basis;
    double multiplicity;
    int maxiter;
    real_t conv;
};

struct rhf_results {
    double energy;
    Matrix F;
    Matrix C;
};

struct uhf_results {
    double energy;
    int nalpha, nbeta;
    Matrix Fa, Fb;
    Matrix Ca, Cb;
};

// Functions
std::vector<libint2::Atom> read_geometry(const std::string &filename);
void print_geometry(const std::vector<libint2::Atom> &atoms);
params read_config(const std::string& config_file);
size_t nbasis(const std::vector<libint2::Shell> &shells);

// Methods
rhf_results RHF(const std::vector<libint2::Atom>& atoms, const libint2::BasisSet& obs, real_t nao, real_t nelectron, params config);
uhf_results UHF(const std::vector<libint2::Atom>& atoms, const libint2::BasisSet& obs, real_t nao, real_t nelectron, params config);

int main(int argc, char *argv[]) {

    using std::cerr;
    using std::cout;
    using std::endl;

    using libint2::BasisSet;

    cout << std::setprecision(12);

    // Reading geometry and config from input files
    std::vector<libint2::Atom> atoms = read_geometry(argv[1]);
    auto config = read_config(argv[2]);

    // Counting the number of electrons
    auto nelectron = 0;
    for (auto i = 0; i < atoms.size(); ++i)
        nelectron += atoms[i].atomic_number;
    const auto ndocc = nelectron / 2;

    print_geometry(atoms);
    cout << "Number of electrons = " << nelectron << endl;


    BasisSet obs(config.basis, atoms);
    auto nao = nbasis(obs.shells());
    cout << endl
         << "Method: " << config.type << endl
         << "Basis Set: " << config.basis << endl
         << "Number of basis functions = " << nao << endl;

    // Main
    if(config.type == "RHF" || config.type == "rhf")
        auto hf_results = RHF(atoms, obs, nao, nelectron, config);
    else if(config.type == "UHF" || config.type == "uhf")
        auto uhf_results = UHF(atoms, obs, nao, nelectron, config);
    else if(config.type == "MP2" || config.type == "mp2")
        ;
    else{
        cout << endl
             << config.type << " is not a supported method" << endl;
    }
}

