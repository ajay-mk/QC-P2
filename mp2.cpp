//
// Created by Ajay Melekamburath on 11/14/22.
//
// Contains functions relevant to MP2 Calculations

#include "mp2.h"


// Function Definitions

DTensor make_so_moes(const Vector &eps_a, const Vector &eps_b, const int &nao) {
    auto n = nao * 2;
    DTensor eps_so(n, n);
    for (auto i = 0; i < n; ++i) {
        if (i % 2 == 0)
            eps_so(i, i) = eps_a(i / 2);
        else if (i % 2 == 1)
            eps_so(i, i) = eps_b(i / 2);
    }
    return eps_so;
}

mp2_output run_mp2(const DTensor &oovv, const DTensor &denom) {
    auto no = oovv.extent(0);
    auto nv = oovv.extent(2);
    mp2_output output;
    real_t energy = 0.0;
    DTensor T(no, no, nv, nv);
    for (auto i = 0; i < no; ++i) {
        for (auto j = 0; j < no; ++j) {
            for (auto a = 0; a < nv; ++a) {
                for (auto b = 0; b < nv; ++b) {
                    T(i, j, a, b) = oovv(i, j, a, b) / denom(i, j, a, b);
                    energy += 0.25 * (oovv(i, j, a, b) * oovv(i, j, a, b)) / denom(i, j, a, b);
                }
            }
        }
    }
    output.T = T;
    output.energy = energy;
    return output;
}

// Main MP2 Function
mp2_results MP2(const libint2::BasisSet &obs, const scf_results &scf, const params &config) {
    mp2_results results;
    std::cout << std::endl
              << "Starting MP2 calculation" << std::endl
              << std::endl;
    if (config.ref == "RHF" || config.ref == "rhf") {
        // Calculating AO Integrals
        auto ao_ints = eri_ao_tensor(obs);

        // Transform AO to spatial MO basis
        auto mo_ints = transform_ao_mo(ao_ints, scf.C, scf.C);

        // Transform to spin MO basis
        results.so_ints = transform_to_so(mo_ints, mo_ints, mo_ints);

        auto moes = make_so_moes(scf.moes, scf.moes, scf.nao);

        auto oovv = slice_ints(results.so_ints, scf.no, scf.nv, "oovv");//(ij|ab)

        auto denom = make_denom(moes, scf.no, scf.nv);

        results.energy = run_mp2(oovv, denom).energy;
        results.T = run_mp2(oovv, denom).T;
    }

    else if (config.ref == "UHF" || config.ref == "uhf") {
        auto ao_ints = eri_ao_tensor(obs);

        auto mo_ints_aa = transform_ao_mo(ao_ints, scf.Ca, scf.Ca);
        auto mo_ints_bb = transform_ao_mo(ao_ints, scf.Cb, scf.Cb);
        auto mo_ints_ab = transform_ao_mo(ao_ints, scf.Ca, scf.Cb);

        results.so_ints = transform_to_so(mo_ints_aa, mo_ints_bb, mo_ints_ab);
        auto oovv = slice_ints(results.so_ints, scf.no, scf.nv, "oovv");//(ij|ab)

        auto moes = make_so_moes(scf.moes_a, scf.moes_b, scf.nao);
        auto denom = make_denom(moes, scf.no, scf.nv);

        results.energy = run_mp2(oovv, denom).energy;
        results.T = run_mp2(oovv, denom).T;
    }
    std::cout << "MP2 Energy: " << results.energy << " Eh" << std::endl;
    return results;
}

// EOF
