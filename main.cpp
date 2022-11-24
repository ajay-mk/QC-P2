#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigenvalues>

// Libint Gaussian integrals library
#include <libint2.hpp>
#if !LIBINT2_CONSTEXPR_STATICS
#include <libint2/statics_definition.h>
#endif
#include "btas/btas.h"
// Include Headers
#include "general.h"
#include "hf.h"

// Typedefs
typedef btas::Tensor<double> DTensor;

struct mp2_results{
    real_t energy;
    DTensor T;
};

// Functions
size_t nbasis(const std::vector<libint2::Shell> &shells);

// Methods
mp2_results MP2(const libint2::BasisSet& obs, const scf_results& scf);
int main(int argc, char *argv[]) {

    using std::cerr;
    using std::cout;
    using std::endl;

    using libint2::BasisSet;

    cout << std::setprecision(12);

    // Reading geometry and config from input files
    auto config = read_config(argv[1]);
    std::vector<libint2::Atom> atoms = read_geometry(config.inputfile);

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
        auto hf_result = RHF(atoms, obs, nao, nelectron, config);
    else if(config.type == "UHF" || config.type == "uhf")
        auto uhf_result = UHF(atoms, obs, nao, nelectron, config);
    else if((config.type == "MP2" || config.type == "mp2") && (config.scf == "RHF" || config.scf == "rhf")){
        auto hf_result = RHF(atoms, obs, nao, nelectron, config);
        auto mp2_result = MP2(obs, hf_result);
        cout << "MP2 Corrected energy: " << hf_result.energy + mp2_result.energy << " Eh" << endl;
    }
    else if((config.type == "MP2" || config.type == "mp2") && (config.scf == "UHF" || config.scf == "UHF")){
        cout << "UMP2 Under development" << endl;
        //auto hf_result = UHF(atoms, obs, nao, nelectron, config);
        //auto mp2_result = MP2(obs, hf_result);
        //cout << "MP2 Corrected energy: " << hf_result.energy + mp2_result.energy << " Eh" << endl;
    }
    else{
        cout << endl
             << config.type << " is not a supported method" << endl;
    }
}

// EOF