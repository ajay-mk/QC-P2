//
// Created by Ajay Melekamburath on 11/14/22.
//
// Contains functions relevant to MP2 Calculations

#include "mp2.h"


// Function Definitions

Vector make_so_moes(const Vector &eps, const int& nao) {
    auto n = nao * 2;
    Vector eps_so(n);
    for (auto i = 0; i < n; i++) {
        eps_so(i) = eps(i / 2);
    }
    return eps_so;
}

Vector make_so_moes_uhf(const Vector& eps_a, const Vector& eps_b, const int& nao){
    auto n = nao * 2;
    Vector eps_so(n);
    for (auto i = 0; i < n; i++) {
        if (i%2 == 0)
            eps_so(i) = eps_a(i/2);
        else if (i%2 == 1)
            eps_so(i) = eps_b(i/2);
    }
    return eps_so;
}


DTensor get_ijab(const DTensor &so_ints, const int &noo, const int &nvo) {
    DTensor ij_ab(noo, noo, nvo, nvo);
    for (auto i = 0; i < noo; i++) {
        for (auto j = 0; j < noo; j++) {
            for (auto a = 0; a < nvo; a++) {
                for (auto b = 0; b < nvo; b++) {
                    ij_ab(i, j, a, b) = so_ints(i, j, noo + a, noo + b);
                }
            }
        }
    }
    return ij_ab;
}

real_t mp2_energy(DTensor ij_ab, Vector eps_so) {
    real_t energy;
    auto noo = ij_ab.extent(0);
    auto nvo = ij_ab.extent(2);
    energy = 0.0;
    for (auto i = 0; i < noo; i++) {
        for (auto j = 0; j < noo; j++) {
            for (auto a = 0; a < nvo; a++) {
                for (auto b = 0; b < nvo; b++) {
                    energy += 0.25 * ij_ab(i, j, a, b) * ij_ab(i, j, a, b) / (eps_so[i] + eps_so[j] - eps_so[noo + a] - eps_so[noo + b]);
                }
            }
        }
    }
    return energy;
}


// Main MP2 Function
mp2_results MP2(const libint2::BasisSet &obs, const scf_results &scf, const params &config) {
    mp2_results results;
    std::cout << std::endl
              << "Starting MP2 calculation" << std::endl
              << std::endl;
    if (config.scf == "RHF" || config.scf == "rhf") {
        // Calculating AO Integrals
        auto ao_ints = eri_ao_tensor(obs);
        // Transform AO to spatial MO basis
        auto mo_ints = transform_ao_mo(ao_ints, scf.C);
        // Transform to spin MO basis
        auto so_ints = transform_to_so(mo_ints);
        auto moes = make_so_moes(scf.moes, scf.nao);
        auto ij_ab = get_ijab(so_ints, scf.noo, scf.nvo);
        results.energy = mp2_energy(ij_ab, moes);
        results.T = ij_ab;
    }
    else if (config.scf == "UHF" || config.scf == "uhf") {
        auto ao_ints = eri_ao_tensor(obs);
        auto mo_ints_aa = transform_ao_mo_uhf(ao_ints, scf.Ca, scf.Ca);
        auto mo_ints_bb = transform_ao_mo_uhf(ao_ints, scf.Cb, scf.Cb);
        auto mo_ints_ab = transform_ao_mo_uhf(ao_ints, scf.Ca, scf.Cb);
        auto so_ints = transform_to_so_uhf(mo_ints_aa, mo_ints_bb, mo_ints_ab);
        auto moes = make_so_moes_uhf(scf.moes_a, scf.moes_b, scf.nao);
        //std::cout << moes << std::endl;
        auto ij_ab = get_ijab(so_ints, scf.noo, scf.nvo);
        results.energy = mp2_energy(ij_ab, moes);
        results.T = ij_ab;
    }
    std::cout << "MP2 Energy: " << results.energy << " Eh" << std::endl;
    return results;
}
// EOF