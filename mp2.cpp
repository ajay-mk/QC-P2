//
// Created by Ajay Melekamburath on 11/14/22.
//
// Contains functions relevant to MP2 Calculations

// Libint Gaussian integrals library
#include <libint2.hpp>
#include <libint2/chemistry/sto3g_atomic_density.h>
#if !LIBINT2_CONSTEXPR_STATICS
#  include <libint2/statics_definition.h>
#endif
#include "btas/btas.h"

#include "mp2.h"

// Function Definitions
// Integral engines are from the example: https://github.com/evaleev/libint/blob/master/tests/hartree-fock/hartree-fock.cc
DTensor eri_ao_tensor(const libint2::BasisSet& obs) {
    using libint2::Shell;
    using libint2::Engine;
    using libint2::Operator;

    const auto n = nbasis(obs.shells());
    DTensor ao_ints(n, n, n, n);
    ao_ints.fill(0.0);

    libint2::initialize();

    Engine engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 0);

    auto shell2bf = obs.shell2bf();

    const auto &buf = engine.results();

    for (auto s1 = 0; s1 != obs.shells().size(); ++s1) {

        auto bf1_first = shell2bf[s1];    // first basis function in this shell
        auto n1 = obs.shells()[s1].size();// number of basis functions in this shell

        for (auto s2 = 0; s2 <= s1; ++s2) {

            auto bf2_first = shell2bf[s2];
            auto n2 = obs.shells()[s2].size();

            for (auto s3 = 0; s3 <= s1; ++s3) {

                auto bf3_first = shell2bf[s3];
                auto n3 = obs.shells()[s3].size();

                const auto s4_max = (s1 == s3) ? s2 : s3;
                for (auto s4 = 0; s4 <= s4_max; ++s4) {

                    auto bf4_first = shell2bf[s4];
                    auto n4 = obs.shells()[s4].size();

                    // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
                    auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
                    auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
                    auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
                    auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

                    engine.compute(obs.shells()[s1], obs.shells()[s2], obs.shells()[s3], obs.shells()[s4]);
                    const auto *buf_1234 = buf[0];
                    if (buf_1234 == nullptr)
                        continue;// if all integrals screened out, skip to next quartet

                    for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
                        const auto bf1 = f1 + bf1_first;
                        for (auto f2 = 0; f2 != n2; ++f2) {
                            const auto bf2 = f2 + bf2_first;
                            for (auto f3 = 0; f3 != n3; ++f3) {
                                const auto bf3 = f3 + bf3_first;
                                for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                                    const auto bf4 = f4 + bf4_first;

                                    const auto value = buf_1234[f1234];

                                    ao_ints(bf1, bf2, bf3, bf4) = value;
                                    ao_ints(bf1, bf2, bf4, bf3) = value;
                                    ao_ints(bf2, bf1, bf3, bf4) = value;
                                    ao_ints(bf2, bf1, bf4, bf3) = value;
                                    ao_ints(bf3, bf4, bf2, bf1) = value;
                                    ao_ints(bf3, bf4, bf1, bf2) = value;
                                    ao_ints(bf4, bf3, bf1, bf2) = value;
                                    ao_ints(bf4, bf3, bf2, bf1) = value;

                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return ao_ints;
}

// AO to MO transformation function : (pq|rs) --> (ia|jb) in Chemist's notation
// <pr|qs> --> (ij|ab) in Dirac notation
DTensor transform_ao_mo(const DTensor& pq_rs, const Matrix& C){
    using btas::contract;
    DTensor ia_jb;
    const int n = pq_rs.extent(0);

    DTensor C0(n,n);
    for (auto a = 0; a < n; a++){
        for (auto b = 0; b < n; b++){
            C0(a, b) = C(a, b);
        }
    }
    // Tensor Contractions
    DTensor pq_rl(n, n, n, n), pq_kl(n, n, n, n), pj_kl(n, n, n, n), ij_kl(n, n, n, n);
    enum {p, q, r, s, i, j, k, l};
    // sum{s} C_{s}^{b} (pq|rs)
    contract(1.0, pq_rs, {p, q, r, s}, C0, {s, l}, 1.0, pq_rl, {p, q, r, l});
    // sum{r} C_{r}^{j} (sum{s} C_{s}^{b} (pq|rs))
    contract(1.0, pq_rl, {p, q, r, l}, C0, {r, k}, 1.0, pq_kl, {p, q, k, l});
    // sum{q} C_{q}^{a} (sum{r} C_{r}^{j} (sum{s} C_{s}^{b} (pq|rs)))
    contract(1.0, pq_kl, {p, q, k, l}, C0, {q, j}, 1.0, pj_kl, {p, j, k, l});
    // sum{p} C_{p}^{i} (sum{q} C_{q}^{a} (sum{r} C_{r}^{j} (sum{s} C_{s}^{b} (pq|rs))))
    contract(1.0, pj_kl, {p, j, k, l}, C0, {p, i}, 1.0, ij_kl, {i, j, k, l});

    //don't need other three tensors anymore
    pq_kl(0,0,0,0), pj_kl(0,0,0,0), pq_rl(0,0,0,0);
    return ij_kl;
}

DTensor transform_ao_mo_uhf(const DTensor& pq_rs, const Matrix& Coeff1, const Matrix& Coeff2){
    using btas::contract;
    DTensor ia_jb;
    const int n = pq_rs.extent(0);

    DTensor Ca(n, n);
    DTensor Cb(n, n);
    for (auto a = 0; a < n; a++){
        for (auto b = 0; b < n; b++){
            Ca(a, b) = Coeff1(a, b);
            Cb(a, b) = Coeff2(a, b);
        }
    }
    // Tensor Contractions
    DTensor pq_rl(n, n, n, n), pq_kl(n, n, n, n), pj_kl(n, n, n, n), ij_kl(n, n, n, n);
    enum {p, q, r, s, i, j, k, l};
    // sum{s} C_{s}^{b} (pq|rs)
    contract(1.0, pq_rs, {p, q, r, s}, Ca, {s, l}, 1.0, pq_rl, {p, q, r, l});
    // sum{r} C_{r}^{j} (sum{s} C_{s}^{b} (pq|rs))
    contract(1.0, pq_rl, {p, q, r, l}, Ca, {r, k}, 1.0, pq_kl, {p, q, k, l});
    // sum{q} C_{q}^{a} (sum{r} C_{r}^{j} (sum{s} C_{s}^{b} (pq|rs)))
    contract(1.0, pq_kl, {p, q, k, l}, Cb, {q, j}, 1.0, pj_kl, {p, j, k, l});
    // sum{p} C_{p}^{i} (sum{q} C_{q}^{a} (sum{r} C_{r}^{j} (sum{s} C_{s}^{b} (pq|rs))))
    contract(1.0, pj_kl, {p, j, k, l}, Cb, {p, i}, 1.0, ij_kl, {i, j, k, l});

    //don't need other three tensors anymore
    pq_kl(0,0,0,0), pj_kl(0,0,0,0), pq_rl(0,0,0,0);
    return ij_kl;
}


DTensor transform_to_so(const DTensor& mo_ints){
    auto n = mo_ints.extent(0) * 2;
    DTensor so_ints(n, n, n, n);
    for (auto i = 0; i < n; i++){
        for (auto j = 0; j < n; j++){
            for (auto k = 0; k < n; k++){
                for (auto l = 0; l < n; l++){
                    so_ints(i, j, k, l) = mo_ints(i/2, k/2, j/2, l/2) * (i%2 == k%2) * (j%2 == l%2)
                                          - mo_ints(j/2, k/2, i/2, l/2) * (j%2 == k%2) * (i%2 == l%2);
                }
            }
        }
    }
    return so_ints;
}

DTensor transform_to_so_uhf(const DTensor& mo_ints_aa, const DTensor& mo_ints_bb, const DTensor& mo_ints_ab){
    auto n = mo_ints_aa.extent(0) * 2;
    DTensor so_ints(n, n, n, n);
    for (auto i = 0; i < n; i++) {
        for (auto j = 0; j < n; j++) {
            for (auto k = 0; k < n; k++) {
                for (auto l = 0; l < n; l++) {
                    if (i % 2 == 0 && k % 2 == 0 && j % 2 == 0 && l % 2 == 0)
                        so_ints(i, j, k, l) = mo_ints_aa(i / 2, k / 2, j / 2, l / 2) - mo_ints_aa(j / 2, k / 2, i / 2, l / 2);
                    else if (i % 2 == 1 && k % 2 == 1 && j % 2 == 1 && l % 2 == 1)
                        so_ints(i, j, k, l) = mo_ints_bb(i / 2, k / 2, j / 2, l / 2) - mo_ints_bb(j / 2, k / 2, i / 2, l / 2);
                    else if (i % 2 == 0 && k % 2 == 0 && j % 2 == 1 && l % 2 == 1)
                        so_ints(i, j, k, l) = mo_ints_ab(i / 2, k / 2, j / 2, l / 2);
                    else if (i % 2 == 1 && k % 2 == 1 && j % 2 == 0 && l % 2 == 0)
                        so_ints(i, j, k, l) = mo_ints_ab(i / 2, k / 2, j / 2, l / 2);
                    else if (i % 2 == 1 && k % 2 == 0 && j % 2 == 0 && l % 2 == 1)
                        so_ints(i, j, k, l) = -mo_ints_ab(j / 2, k / 2, i / 2, l / 2);
                    else if (i % 2 == 0 && k % 2 == 1 && j % 2 == 1 && l % 2 == 0)
                        so_ints(i, j, k, l) = -mo_ints_ab(j / 2, k / 2, i / 2, l / 2);
                }
            }
        }
    }
    return so_ints;
    }

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

//DTensor get_iajb(const DTensor& ij_kl, const int& nocc, const int& nuocc)
//{
//    DTensor ia_jb(nocc, nuocc, nocc, nuocc);
//    auto n = ij_kl.extent(0);
//    for (auto i = 0; i < nocc; i++){
//        for (auto a = 0; a < nuocc; a++){
//            for (auto j = 0; j < nocc; j++){
//                for (auto b = 0; b < nuocc; b++){
//                    ia_jb(i,a,j,b) = ij_kl(i,nocc+a,j,nocc+b);
//                }
//            }
//        }
//    }
//    return ia_jb;
//}

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

// Calculates MP2 Energy - MO Basis
//real_t mp2_energy(DTensor ia_jb, Vector eps){
//    real_t energy;
//    auto nocc = ia_jb.extent(0);
//    auto nuocc = ia_jb.extent(1);
//    energy = 0.0;
//    for (auto i = 0; i < nocc; i++){
//        for (auto a = 0; a < nuocc; a++){
//            for (auto j = 0; j < nocc; j++){
//                for (auto b = 0; b < nuocc; b++){
//                    energy += ia_jb(i, a, j, b) *
//                              (2 * ia_jb(i, a, j, b) - ia_jb(i, b, j, a))/
//                              (eps[i] + eps[j] - eps[nocc + a] - eps[nocc + b]);
//                }
//            }
//        }
//    }
//    return energy;
//}

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
        //results.T = ij_ab;
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
        //results.T = ij_ab;
    }
    std::cout << "MP2 Energy: " << results.energy << " Eh" << std::endl;
    return results;
}
// EOF