#include <iomanip>
#include <iostream>
#include <vector>


// Include Headers
#include "cc.h"
#include "general.h"
#include "hf.h"
#include "mp2.h"

// Main
int main(int argc, char *argv[]) {


    using std::cerr;
    using std::cout;
    using std::endl;

    // Setting cout precision
    cout << std::setprecision(12);

    using libint2::BasisSet;

    // Reading geometry and config from input files
    auto config = read_config_json(argv[1]);
    std::vector<libint2::Atom> atoms = read_geometry(config.inputfile);

    // Counting the number of electrons
    auto nelectron = 0;
    for (auto &atom: atoms)
        nelectron += atom.atomic_number;

    print_geometry(atoms);
    cout << "Number of electrons = " << nelectron << endl;


    BasisSet obs(config.basis, atoms);
    auto nao = nbasis(obs.shells());
    cout << endl
         << "Method: " << config.type << endl
         << "Reference: " << config.ref << endl
         << "Basis Set: " << config.basis << endl
         << "Number of basis functions = " << nao << endl;

    // Main

    // HF Bracket
    if (config.type == "RHF" || config.type == "rhf")
        auto hf_result = RHF(atoms, obs, nelectron, config);
    else if (config.type == "UHF" || config.type == "uhf")
        auto hf_result = UHF(atoms, obs, nelectron, config);

    // MP2 Bracket
    else if (config.type == "MP2" || config.type == "mp2") {
        scf_results hf_result;
        if (config.ref == "RHF" || config.ref == "rhf")
            hf_result = RHF(atoms, obs, nelectron, config);
        else if (config.ref == "UHF" || config.ref == "UHF")
            hf_result = UHF(atoms, obs, nelectron, config);
        auto mp2_result = MP2(obs, hf_result, config);

        cout << "Total MP2 energy: " << hf_result.energy + mp2_result.energy << " Eh" << endl;
    }

    // Coupled Cluster Bracket
    else if (config.type == "CCSD" || config.type == "ccsd") {
        scf_results hf_result;
        if (config.ref == "RHF" || config.ref == "rhf")
            hf_result = RHF(atoms, obs, nelectron, config);
        else if (config.ref == "UHF" || config.ref == "UHF")
            hf_result = UHF(atoms, obs, nelectron, config);
        auto mp2_result = MP2(obs, hf_result, config);
        auto ccsd_result = CCSD(hf_result, mp2_result, config);

        cout << "Total energy: " << hf_result.energy + ccsd_result.ccsd_energy << " Eh" << endl;
    }

    else if (config.type == "CCSD(T)" || config.type == "ccsd(t)") {
        scf_results hf_result;
        if (config.ref == "RHF" || config.ref == "rhf")
            hf_result = RHF(atoms, obs, nelectron, config);
        else if (config.ref == "UHF" || config.ref == "UHF")
            hf_result = UHF(atoms, obs, nelectron, config);
        auto mp2_result = MP2(obs, hf_result, config);
        auto ccsd_result = CCSD(hf_result, mp2_result, config);
        auto moes = make_moe_tensors(hf_result, config);
        auto t_energy = CCSD_T(ccsd_result, moes);

        cout << "Total energy: " << hf_result.energy + ccsd_result.ccsd_energy + t_energy << " Eh" << endl;

    }
    // Other Methods
    else {
        cout << endl
             << config.type << " is not a supported method" << endl;
    }
}
// EOF
