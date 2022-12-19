#include <iostream>
#include <vector>

// Libint Gaussian integrals library
#include <libint2.hpp>
#if !LIBINT2_CONSTEXPR_STATICS
#include <libint2/statics_definition.h>
#endif

// Include Headers
#include "general.h"
#include "hf.h"
#include "mp2.h"
#include "cc.h"

// Main
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
    for (auto & atom : atoms)
        nelectron += atom.atomic_number;

    print_geometry(atoms);
    cout << "Number of electrons = " << nelectron << endl;


    BasisSet obs(config.basis, atoms);
    auto nao = nbasis(obs.shells());
    cout << endl
         << "Method: " << config.type << endl
         << "SCF Type: " << config.scf << endl
         << "Basis Set: " << config.basis << endl
         << "Number of basis functions = " << nao << endl;

    // Main

    // HF Bracket
    if(config.type == "RHF" || config.type == "rhf")
        auto hf_result = RHF(atoms, obs, nao, nelectron, config);
    else if(config.type == "UHF" || config.type == "uhf")
        auto uhf_result = UHF(atoms, obs, nao, nelectron, config);
    // MP2 Bracket
    else if((config.type == "MP2" || config.type == "mp2") && (config.scf == "RHF" || config.scf == "rhf")){
        auto hf_result = RHF(atoms, obs, nao, nelectron, config);
        auto mp2_result = MP2(obs, hf_result, config);
        cout << "Total MP2 energy: " << hf_result.energy + mp2_result.energy << " Eh" << endl;
    }
    else if((config.type == "MP2" || config.type == "mp2") && (config.scf == "UHF" || config.scf == "UHF")){
        auto hf_result = UHF(atoms, obs, nao, nelectron, config);
        auto mp2_result = MP2(obs, hf_result, config);
        cout << "Total MP2 energy: " << hf_result.energy + mp2_result.energy << " Eh" << endl;
    }
    // CCSD Bracket
    else if((config.type == "CCSD" || config.type == "ccsd") && (config.scf == "RHF" || config.scf == "rhf")){
        auto hf_result = RHF(atoms, obs, nao, nelectron, config);
        auto mp2_result = MP2(obs, hf_result, config);
        auto CC_Energy = CCSD(hf_result, mp2_result, config);
    }
    else if((config.type == "CCSD" || config.type == "ccsd") && (config.scf == "UHF" || config.scf == "UHF")) {
        auto hf_result = UHF(atoms, obs, nao, nelectron, config);
        auto mp2_result = MP2(obs, hf_result, config);
        auto CC_Energy = CCSD(hf_result, mp2_result, config);
    }
    // Other Methods
    else{
        cout << endl
             << config.type << " is not a supported method" << endl;
    }
}
// EOF