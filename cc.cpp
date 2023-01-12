//
// Created by Ajay Melekamburath on 12/5/22.
//
// // Contains functions relevant to CCSD & CCSD(T)

#include "cc.h"

int_struct get_integrals(const scf_results &SCF, const mp2_results &MP2) {
    int_struct output;
    output.oooo = slice_ints(MP2.so_ints, SCF.no, SCF.nv, "oooo");
    output.ovoo = slice_ints(MP2.so_ints, SCF.no, SCF.nv, "ovoo");
    output.oovo = slice_ints(MP2.so_ints, SCF.no, SCF.nv, "oovo");
    output.ooov = slice_ints(MP2.so_ints, SCF.no, SCF.nv, "ooov");
    output.ovov = slice_ints(MP2.so_ints, SCF.no, SCF.nv, "ovov");
    output.ovvo = slice_ints(MP2.so_ints, SCF.no, SCF.nv, "ovvo");
    output.oovv = slice_ints(MP2.so_ints, SCF.no, SCF.nv, "oovv");
    output.vovv = slice_ints(MP2.so_ints, SCF.no, SCF.nv, "vovv");
    output.ovvv = slice_ints(MP2.so_ints, SCF.no, SCF.nv, "ovvv");
    output.vvvo = slice_ints(MP2.so_ints, SCF.no, SCF.nv, "vvvo");
    output.vvvv = slice_ints(MP2.so_ints, SCF.no, SCF.nv, "vvvv");

    return output;
}

DTensor make_fock(const Vector &eps1, const Vector &eps2, int n1, int n2) {
    DTensor result(n2 - n1, n2 - n1);
    for (auto i = 0; i < n2 - n1; i++) {
        if (n1 % 2 == 0)
            result(i, i) = (i % 2 == 0) * eps1(n1 / 2 + i / 2) + (i % 2 == 1) * eps2(n1 / 2 + i / 2);
        else if (n1 % 2 == 1)
            result(i, i) = (i % 2 == 0) * eps1(n1 / 2 + i / 2 + 1) + (i % 2 == 1) * eps2(n1 / 2 + i / 2 + 1);
    }
    return result;
}

moes make_moe_tensors(const scf_results &scf, const params &config) {
    moes result;

    if (config.ref == "RHF") {
        result.F = make_fock(scf.moes, scf.moes, 0, scf.no + scf.nv);
        result.F_ii = make_fock(scf.moes, scf.moes, 0, scf.no);
        result.F_aa = make_fock(scf.moes, scf.moes, scf.no, scf.no + scf.nv);
        DTensor F_ia(scf.no, scf.nv);
        for (auto i = 0; i < scf.no; i++) {
            for (auto a = 0; a < scf.nv; a++) {
                F_ia(i, a) = result.F(i, scf.no + a);
            }
        }
        result.F_ia = F_ia;
    }

    if (config.ref == "UHF") {
        result.F = make_fock(scf.moes_a, scf.moes_b, 0, scf.no + scf.nv);
        result.F_ii = make_fock(scf.moes_a, scf.moes_b, 0, scf.no);

        DTensor F_aa(scf.nv, scf.nv);
        for (auto a = 0; a < scf.nv; a++) {
            F_aa(a, a) = result.F(scf.no + a, scf.no + a);
        }
        result.F_aa = F_aa;

        DTensor F_ia(scf.no, scf.nv);
        for (auto i = 0; i < scf.no; i++) {
            for (auto a = 0; a < scf.nv; a++) {
                F_ia(i, a) = result.F(i, scf.no + a);
            }
        }
    }
    return result;
}

// Stanton Equation #12
DTensor make_D_ia(const moes &moes) {
    auto no = moes.F_ii.extent(0);
    auto nv = moes.F_aa.extent(0);
    DTensor D_ia(no, nv);
    for (auto i = 0; i < no; i++) {
        for (auto a = 0; a < nv; a++) {
            D_ia(i, a) = moes.F_ii(i, i) - moes.F_aa(a, a);
        }
    }
    return D_ia;
}

// Stanton Equation #13
DTensor make_D_ijab(const moes &moes) {
    auto no = moes.F_ii.extent(0);
    auto nv = moes.F_aa.extent(0);
    DTensor D_ijab(no, no, nv, nv);
    for (auto i = 0; i < no; i++) {
        for (auto j = 0; j < no; j++) {
            for (auto a = 0; a < nv; a++) {
                for (auto b = 0; b < nv; b++) {
                    D_ijab(i, j, a, b) = moes.F_ii(i, i) + moes.F_ii(j, j) - moes.F_aa(a, a) - moes.F_aa(b, b);
                }
            }
        }
    }
    return D_ijab;
}

DTensor make_D_triples(const moes &moes) {
    auto no = moes.F_ii.extent(0);
    auto nv = moes.F_aa.extent(1);
    DTensor D_triples(no, no, no, nv, nv, nv);
    D_triples.fill(0.0);
    for (auto i = 0; i < no; i++) {
        for (auto j = 0; j < no; j++) {
            for (auto k = 0; k < no; k++) {
                for (auto a = 0; a < nv; a++) {
                    for (auto b = 0; b < nv; b++) {
                        for (auto c = 0; c < nv; c++) {
                            D_triples(i, j, k, a, b, c) = moes.F_ii(i, i) + moes.F_ii(j, j) + moes.F_ii(k, k) - moes.F_aa(a, a) - moes.F_aa(b, b) - moes.F_aa(c, c);
                        }
                    }
                }
            }
        }
    }
    return D_triples;
}

DTensor make_tau(const DTensor &Ts, const DTensor &Td) {
    auto no = Ts.extent(0);
    auto nv = Ts.extent(1);
    DTensor tau(no, no, nv, nv);
    for (auto i = 0; i < no; i++) {
        for (auto j = 0; j < no; j++) {
            for (auto a = 0; a < nv; a++) {
                for (auto b = 0; b < nv; b++) {
                    tau(i, j, a, b) = Td(i, j, a, b) + Ts(i, a) * Ts(j, b) - Ts(i, b) * Ts(j, a);// Stanton Equation #9
                }
            }
        }
    }
    return tau;
}
DTensor make_tau_bar(const DTensor &Ts, const DTensor &Td) {
    auto no = Ts.extent(0);
    auto nv = Ts.extent(1);
    DTensor tau_bar(no, no, nv, nv);
    for (auto i = 0; i < no; i++) {
        for (auto j = 0; j < no; j++) {
            for (auto a = 0; a < nv; a++) {
                for (auto b = 0; b < nv; b++) {
                    tau_bar(i, j, a, b) = Td(i, j, a, b) +
                                          0.5 * (Ts(i, a) * Ts(j, b) - Ts(i, b) * Ts(j, a));// Stanton Equation #10
                }
            }
        }
    }
    return tau_bar;
}

DTensor multiply_Ts(const DTensor &Ts) {
    auto no = Ts.extent(0);
    auto nv = Ts.extent(1);
    DTensor result(no, nv, no, nv);
    for (auto j = 0; j < no; j++) {
        for (auto n = 0; n < no; n++) {
            for (auto f = 0; f < nv; f++) {
                for (auto b = 0; b < nv; b++) {
                    result(j, f, n, b) = Ts(j, f) * Ts(n, b);
                }
            }
        }
    }
    return result;
}

// Intermediates Update
cc_intermediates update_intermediates(const DTensor &Ts, const DTensor &Td,
                                      const int_struct &integrals, const moes &moes) {
    using btas::contract;
    cc_intermediates intermediates;
    auto no = Ts.extent(0);
    auto nv = Ts.extent(1);
    enum { a,
           b,
           i,
           j,
           m,
           n,
           e,
           f };

    DTensor F_ae(nv, nv);// Stanton Equation #3

    contract(-0.5, moes.F_ia, {m, e}, Ts, {m, a}, 1.0, F_ae, {a, e});
    contract(1.0, Ts, {m, f}, integrals.ovvv, {m, a, f, e}, 1.0, F_ae, {a, e});
    contract(-0.5, make_tau_bar(Ts, Td), {m, n, a, f}, integrals.oovv, {m, n, e, f}, 1.0, F_ae, {a, e});

    for (auto a = 0; a < nv; a++) {
        for (auto e = 0; e < nv; e++) {
            F_ae(a, e) += (1 - (a == e)) * moes.F_aa(a, e);// First term of Equation #3
        }
    }
    intermediates.F_ae = F_ae;// Storing F_ae


    DTensor F_mi(no, no);// Stanton Equation #4
    contract(0.5, Ts, {i, e}, moes.F_ia, {m, e}, 1.0, F_mi, {m, i});
    contract(1.0, Ts, {n, e}, integrals.ooov, {m, n, i, e}, 1.0, F_mi, {m, i});
    contract(0.5, make_tau_bar(Ts, Td), {i, n, e, f}, integrals.oovv, {m, n, e, f}, 1.0, F_mi, {m, i});

    for (auto m = 0; m < no; m++) {
        for (auto i = 0; i < no; i++) {
            F_mi(m, i) += (1 - (m == i)) * moes.F_ii(m, i);// First term of Equation #4
        }
    }
    intermediates.F_mi = F_mi;

    DTensor F_me(no, nv);// Stanton Equation #5
    contract(1.0, Ts, {n, f}, integrals.oovv, {m, n, e, f}, 1.0, F_me, {m, e});
    for (auto m = 0; m < no; m++) {
        for (auto e = 0; e < nv; e++) {
            F_me(m, e) += moes.F_ia(m, e);// First term of Equation #5
        }
    }
    intermediates.F_me = F_me;

    DTensor W_mnij(no, no, no, no);// Stanton Equation #6
    contract(1.0, Ts, {j, e}, integrals.ooov, {m, n, i, e}, 1.0, W_mnij, {m, n, i, j});
    contract(-1.0, Ts, {i, e}, integrals.ooov, {m, n, j, e}, 1.0, W_mnij, {m, n, i, j});
    contract(0.25, make_tau(Ts, Td), {i, j, e, f}, integrals.oovv, {m, n, e, f}, 1.0, W_mnij, {m, n, i, j});
    W_mnij += integrals.oooo;// First term of Equation #6

    intermediates.W_mnij = W_mnij;

    DTensor W_abef(nv, nv, nv, nv);// Stanton Equation #7
    contract(-1.0, Ts, {m, b}, integrals.vovv, {a, m, e, f}, 1.0, W_abef, {a, b, e, f});
    contract(1.0, Ts, {m, a}, integrals.vovv, {b, m, e, f}, 1.0, W_abef, {a, b, e, f});
    contract(0.25, make_tau(Ts, Td), {m, n, a, b}, integrals.oovv, {m, n, e, f}, 1.0, W_abef, {a, b, e, f});
    W_abef += integrals.vvvv;// First term of Equation #7

    intermediates.W_abef = W_abef;

    DTensor W_mbej(no, nv, nv, no);// Stanton Equation #8
    contract(1.0, Ts, {j, f}, integrals.ovvv, {m, b, e, f}, 1.0, W_mbej, {m, b, e, j});
    contract(-1.0, Ts, {n, b}, integrals.oovo, {m, n, e, j}, 1.0, W_mbej, {m, b, e, j});
    contract(-0.5, Td, {j, n, f, b}, integrals.oovv, {m, n, e, f}, 1.0, W_mbej, {m, b, e, j});// Split last term into two
    contract(-1.0, multiply_Ts(Ts), {j, f, n, b}, integrals.oovv, {m, n, e, f}, 1.0, W_mbej, {m, b, e, j});
    W_mbej += integrals.ovvo;// First term of Equation #8

    intermediates.W_mbej = W_mbej;

    return intermediates;
}

DTensor make_T1(const DTensor &Ts, const DTensor &Td, const int_struct &integrals,
                const cc_intermediates &intermediates, const DTensor &D_ia, const moes &moes) {
    using btas::contract;
    auto no = Ts.extent(0);
    auto nv = Ts.extent(1);
    enum { a,
           i,
           m,
           n,
           e,
           f };
    DTensor tempT1(no, nv);
    // Stanton Equation #1
    contract(1.0, Ts, {i, e}, intermediates.F_ae, {a, e}, 1.0, tempT1, {i, a});
    contract(-1.0, Ts, {m, a}, intermediates.F_mi, {m, i}, 1.0, tempT1, {i, a});
    contract(1.0, Td, {i, m, a, e}, intermediates.F_me, {m, e}, 1.0, tempT1, {i, a});
    contract(-1.0, Ts, {n, f}, integrals.ovov, {n, a, i, f}, 1.0, tempT1, {i, a});
    contract(-0.5, Td, {i, m, e, f}, integrals.ovvv, {m, a, e, f}, 1.0, tempT1, {i, a});
    contract(-0.5, Td, {m, n, a, e}, integrals.oovo, {n, m, e, i}, 1.0, tempT1, {i, a});

    // For first part of equation #1
    tempT1 += moes.F_ia;

    // Since Equation #1 is for t_ia * D_ia
    DTensor newT1(no, nv);
    for (auto i = 0; i < no; i++) {
        for (auto a = 0; a < nv; a++) {
            newT1(i, a) = tempT1(i, a) / D_ia(i, a);
        }
    }
    tempT1(0, 0);// Free up

    return newT1;
}

DTensor make_T2(const DTensor &Ts, const DTensor &Td, const int_struct &integrals,
                const cc_intermediates &intermediates, const DTensor &D_ijab) {
    using btas::contract;
    auto no = Ts.extent(0);
    auto nv = Ts.extent(1);

    DTensor tempT2(no, no, nv, nv);
    enum { a,
           b,
           i,
           j,
           m,
           n,
           e,
           f };
    // Stanton Equation #2
    tempT2 += integrals.oovv;// <ij||ab>


    // P_(ab) t_ijae F_be
    contract(1.0, Td, {i, j, a, e}, intermediates.F_ae, {b, e}, 1.0, tempT2, {i, j, a, b});
    contract(-1.0, Td, {i, j, b, e}, intermediates.F_ae, {a, e}, 1.0, tempT2, {i, j, a, b});

    // P_(ab) t_ijae 0.5 * t_mb F_me
    // Do the inner contraction first
    DTensor d1(nv, nv);// Temp tensor
    contract(1.0, Ts, {m, b}, intermediates.F_me, {m, e}, 1.0, d1, {b, e});
    // Contraction for rest of the equation
    contract(-0.5, Td, {i, j, a, e}, d1, {b, e}, 1.0, tempT2, {i, j, a, b});
    // Do the same with switching a & b
    d1 -= d1;// Clearing d1
    contract(1.0, Ts, {m, a}, intermediates.F_me, {m, e}, 1.0, d1, {a, e});
    contract(0.5, Td, {i, j, b, e}, d1, {a, e}, 1.0, tempT2, {i, j, a, b});
    d1(0, 0);// Free up


    // -P_(ij) t_imab F_be
    contract(-1.0, Td, {i, m, a, b}, intermediates.F_mi, {m, j}, 1.0, tempT2, {i, j, a, b});
    contract(1.0, Td, {j, m, a, b}, intermediates.F_mi, {m, i}, 1.0, tempT2, {i, j, a, b});


    // -P_(ij) t_imab * t_je F_me
    DTensor d2(no, no);
    contract(1.0, Ts, {j, e}, intermediates.F_me, {m, e}, 1.0, d2, {j, m});
    contract(-0.5, Td, {i, m, a, b}, d2, {j, m}, 1.0, tempT2, {i, j, a, b});
    d2 -= d2;
    contract(1.0, Ts, {i, e}, intermediates.F_me, {m, e}, 1.0, d2, {i, m});
    contract(0.5, Td, {j, m, a, b}, d2, {i, m}, 1.0, tempT2, {i, j, a, b});


    contract(0.5, make_tau(Ts, Td), {m, n, a, b}, intermediates.W_mnij, {m, n, i, j}, 1.0, tempT2, {i, j, a, b});//4th term
    contract(0.5, make_tau(Ts, Td), {i, j, e, f}, intermediates.W_abef, {a, b, e, f}, 1.0, tempT2, {i, j, a, b});//5th term


    //6th term
    // First Part
    contract(1.0, Td, {i, m, a, e}, intermediates.W_mbej, {m, b, e, j}, 1.0, tempT2, {i, j, a, b});
    contract(-1.0, Td, {j, m, a, e}, intermediates.W_mbej, {m, b, e, i}, 1.0, tempT2, {i, j, a, b});
    contract(-1.0, Td, {i, m, b, e}, intermediates.W_mbej, {m, a, e, j}, 1.0, tempT2, {i, j, a, b});
    contract(1.0, Td, {j, m, b, e}, intermediates.W_mbej, {m, a, e, i}, 1.0, tempT2, {i, j, a, b});

    // Second Part
    contract(1.0, multiply_Ts(Ts), {i, e, m, a}, integrals.ovov, {m, b, j, e}, 1.0, tempT2, {i, j, a, b});
    contract(-1.0, multiply_Ts(Ts), {j, e, m, a}, integrals.ovov, {m, b, i, e}, 1.0, tempT2, {i, j, a, b});
    contract(-1.0, multiply_Ts(Ts), {i, e, m, b}, integrals.ovov, {m, a, j, e}, 1.0, tempT2, {i, j, a, b});
    contract(1.0, multiply_Ts(Ts), {j, e, m, b}, integrals.ovov, {m, a, i, e}, 1.0, tempT2, {i, j, a, b});

    //7th term
    contract(1.0, Ts, {i, e}, integrals.vvvo, {a, b, e, j}, 1.0, tempT2, {i, j, a, b});
    contract(-1.0, Ts, {j, e}, integrals.vvvo, {a, b, e, i}, 1.0, tempT2, {i, j, a, b});

    //8th term
    contract(-1.0, Ts, {m, a}, integrals.ovoo, {m, b, i, j}, 1.0, tempT2, {i, j, a, b});
    contract(1.0, Ts, {m, b}, integrals.ovoo, {m, a, i, j}, 1.0, tempT2, {i, j, a, b});

    // Since Equation #2 is for t_ia * D_ijab
    DTensor newT2(no, no, nv, nv);
    for (auto i = 0; i < no; i++) {
        for (auto j = 0; j < no; j++) {
            for (auto a = 0; a < nv; a++) {
                for (auto b = 0; b < nv; b++) {
                    newT2(i, j, a, b) = tempT2(i, j, a, b) / D_ijab(i, j, a, b);
                }
            }
        }
    }
    tempT2(0, 0, 0, 0);// Free up

    return newT2;
}

// TDC & HFS Equation #134
real_t ccsd_energy(const DTensor &ts, const DTensor &td, const DTensor &oovv, const moes &moes) {
    real_t energy = 0;
    auto no = ts.extent(0);
    auto nv = ts.extent(1);
    for (auto i = 0; i < no; i++) {
        for (auto a = 0; a < nv; a++) {
            energy += moes.F_ia(i, a) * ts(i, a);
            for (auto j = 0; j < no; j++) {
                for (auto b = 0; b < nv; b++) {
                    energy += 0.25 * oovv(i, j, a, b) * td(i, j, a, b) + 0.5 * oovv(i, j, a, b) * ts(i, a) * ts(j, b);
                }
            }
        }
    }
    return energy;
}

// CCSD Function

cc_results CCSD(const scf_results &scf, const mp2_results &mp2, const params &config) {
    cc_results results;
    results.sliced_ints = get_integrals(scf, mp2);

    auto moes = make_moe_tensors(scf, config);
    auto D_ia = make_D_ia(moes);
    auto D_ijab = make_D_ijab(moes);

    DTensor Ts(scf.no, scf.nv);
    Ts.fill(0.0);
    results.T1 = Ts;

    results.T2 = mp2.T;// Initial T2 guess from MP2

    // Starting iterations
    real_t E_CC_last = 0.0;
    std::cout << "Starting CCSD Calculations" << std::endl;
    for (auto iter = 0; iter < config.maxiter; iter++) {
        if (iter == 0) {
            std::cout << std::endl
                      << "Iter         E_CC (Eh)        Delta(E_CC)" << std::endl;
        }

        results.ccsd_energy = ccsd_energy(results.T1, results.T2, results.sliced_ints.oovv, moes);
        auto Del_E_CC = results.ccsd_energy - E_CC_last;

        E_CC_last = results.ccsd_energy;

        printf(" %02d %20.12f %20.12f\n", iter, results.ccsd_energy, Del_E_CC);
        if (abs(Del_E_CC) < config.cc_conv) {
            std::cout << "CC energy converged" << std::endl;
            break;
        }
        cc_intermediates intermediates = update_intermediates(results.T1, results.T2, results.sliced_ints, moes);

        // Make T1 & T2
        results.T1 = make_T1(results.T1, results.T2, results.sliced_ints, intermediates, D_ia, moes);
        results.T2 = make_T2(results.T1, results.T2, results.sliced_ints, intermediates, D_ijab);
    }
    std::cout << std::endl
              << "CCSD energy: " << results.ccsd_energy << " Eh" << std::endl;
    return results;
}
// Equations from: https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2306
real_t CCSD_T(const cc_results &ccResults, const moes &moes) {
    using btas::contract;
    auto integrals = ccResults.sliced_ints;
    auto no = ccResults.T1.extent(0);
    auto nv = ccResults.T1.extent(1);

    auto D_triples = make_D_triples(moes);
    enum { i,
           j,
           k,
           a,
           b,
           c,
           e,
           m };


    // Disconnected triples
    //    DTensor tempTd(no, no, no, nv, nv, nv);
    //    contract(1.0, ccResults.T1, {i, a}, integrals.oovv, {j, k, b, c}, 1.0, tempTd, {i, j, k, a, b, c});
    //    contract(1.0, ccResults.T1, {i, b}, integrals.oovv, {j, k, a, c}, -1.0, tempTd, {i, j, k, a, b, c});
    //    contract(1.0, ccResults.T1, {i, c}, integrals.oovv, {j, k, b, a}, -1.0, tempTd, {i, j, k, a, b, c});
    //
    //    contract(1.0, ccResults.T1, {j, a}, integrals.oovv, {i, k, b, c}, -1.0, tempTd, {i, j, k, a, b, c});
    //    contract(1.0, ccResults.T1, {j, b}, integrals.oovv, {i, k, a, c}, 1.0, tempTd, {i, j, k, a, b, c});
    //    contract(1.0, ccResults.T1, {j, c}, integrals.oovv, {i, k, b, a}, 1.0, tempTd, {i, j, k, a, b, c});
    //
    //    contract(1.0, ccResults.T1, {k, a}, integrals.oovv, {j, i, b, c}, -1.0, tempTd, {i, j, k, a, b, c});
    //    contract(1.0, ccResults.T1, {k, b}, integrals.oovv, {j, i, a, c}, -1.0, tempTd, {i, j, k, a, b, c});
    //    contract(1.0, ccResults.T1, {k, c}, integrals.oovv, {j, i, b, a}, -1.0, tempTd, {i, j, k, a, b, c});

    DTensor dT(no, no, no, nv, nv, nv);
    for (auto i = 0; i < no; i++) {
        for (auto j = 0; j < no; j++) {
            for (auto k = 0; k < no; k++) {
                for (auto a = 0; a < nv; a++) {
                    for (auto b = 0; b < nv; b++) {
                        for (auto c = 0; c < nv; c++) {
                            dT(i, j, k, a, b, c) += (ccResults.T1(i, a) * integrals.oovv(j, k, b, c) - ccResults.T1(i, b) * integrals.oovv(j, k, a, c) - ccResults.T1(i, c) * integrals.oovv(j, k, b, a)

                                                     - ccResults.T1(j, a) * integrals.oovv(i, k, b, c) + ccResults.T1(j, b) * integrals.oovv(i, k, a, c) + ccResults.T1(j, c) * integrals.oovv(i, k, b, a)

                                                     - ccResults.T1(k, a) * integrals.oovv(j, i, b, c) + ccResults.T1(k, b) * integrals.oovv(j, i, a, c) + ccResults.T1(k, c) * integrals.oovv(j, i, b, a)) /
                                                    D_triples(i, j, k, a, b, c);

                            //dT(i, j, k, a, b, c) = tempTd(i, j, k, a, b, c)/D_triples(i, j, k, a, b, c);
                        }
                    }
                }
            }
        }
    }
    //tempTd(0, 0, 0, 0, 0, 0);

    // Connected triples - Split into two
    DTensor tempTc(no, no, no, nv, nv, nv);
    contract(1.0, ccResults.T2, {j, k, a, e}, integrals.vovv, {e, i, b, c}, 1.0, tempTc, {i, j, k, a, b, c});
    contract(-1.0, ccResults.T2, {j, k, b, e}, integrals.vovv, {e, i, a, c}, 1.0, tempTc, {i, j, k, a, b, c});
    contract(-1.0, ccResults.T2, {j, k, c, e}, integrals.vovv, {e, i, b, a}, 1.0, tempTc, {i, j, k, a, b, c});

    contract(-1.0, ccResults.T2, {i, k, a, e}, integrals.vovv, {e, j, b, c}, 1.0, tempTc, {i, j, k, a, b, c});
    contract(1.0, ccResults.T2, {i, k, b, e}, integrals.vovv, {e, j, a, c}, 1.0, tempTc, {i, j, k, a, b, c});
    contract(1.0, ccResults.T2, {i, k, c, e}, integrals.vovv, {e, j, b, a}, 1.0, tempTc, {i, j, k, a, b, c});

    contract(-1.0, ccResults.T2, {j, i, a, e}, integrals.vovv, {e, k, b, c}, 1.0, tempTc, {i, j, k, a, b, c});
    contract(1.0, ccResults.T2, {j, i, b, e}, integrals.vovv, {e, k, a, c}, 1.0, tempTc, {i, j, k, a, b, c});
    contract(1.0, ccResults.T2, {j, i, c, e}, integrals.vovv, {e, k, b, a}, 1.0, tempTc, {i, j, k, a, b, c});

    contract(-1.0, ccResults.T2, {i, m, b, c}, integrals.ovoo, {m, a, j, k}, 1.0, tempTc, {i, j, k, a, b, c});
    contract(1.0, ccResults.T2, {i, m, a, c}, integrals.ovoo, {m, b, j, k}, 1.0, tempTc, {i, j, k, a, b, c});
    contract(1.0, ccResults.T2, {i, m, b, a}, integrals.ovoo, {m, c, j, k}, 1.0, tempTc, {i, j, k, a, b, c});

    contract(1.0, ccResults.T2, {j, m, b, c}, integrals.ovoo, {m, a, i, k}, 1.0, tempTc, {i, j, k, a, b, c});
    contract(-1.0, ccResults.T2, {j, m, a, c}, integrals.ovoo, {m, b, i, k}, 1.0, tempTc, {i, j, k, a, b, c});
    contract(-1.0, ccResults.T2, {j, m, b, a}, integrals.ovoo, {m, c, i, k}, 1.0, tempTc, {i, j, k, a, b, c});

    contract(1.0, ccResults.T2, {k, m, b, c}, integrals.ovoo, {m, a, j, i}, 1.0, tempTc, {i, j, k, a, b, c});
    contract(-1.0, ccResults.T2, {k, m, a, c}, integrals.ovoo, {m, b, j, i}, 1.0, tempTc, {i, j, k, a, b, c});
    contract(-1.0, ccResults.T2, {k, m, b, a}, integrals.ovoo, {m, c, j, i}, 1.0, tempTc, {i, j, k, a, b, c});

    DTensor cT(no, no, no, nv, nv, nv);
    for (auto i = 0; i < no; i++) {
        for (auto j = 0; j < no; j++) {
            for (auto k = 0; k < no; k++) {
                for (auto a = 0; a < nv; a++) {
                    for (auto b = 0; b < nv; b++) {
                        for (auto c = 0; c < nv; c++) {
                            cT(i, j, k, a, b, c) = tempTc(i, j, k, a, b, c) / D_triples(i, j, k, a, b, c);
                        }
                    }
                }
            }
        }
    }
    tempTc(0, 0, 0, 0, 0, 0);

    real_t val = 0.0;
    for (auto i = 0; i < no; i++) {
        for (auto j = 0; j < no; j++) {
            for (auto k = 0; k < no; k++) {
                for (auto a = 0; a < nv; a++) {
                    for (auto b = 0; b < nv; b++) {
                        for (auto c = 0; c < nv; c++) {
                            val += (cT(i, j, k, a, b, c) * D_triples(i, j, k, a, b, c)) * (cT(i, j, k, a, b, c) + dT(i, j, k, a, b, c));
                        }
                    }
                }
            }
        }
    }
    auto t_energy = val / 36;

    std::cout << "(T) energy: " << t_energy << " Eh" << std::endl;
    return t_energy;
}

// EOF