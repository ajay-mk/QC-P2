//
// Created by Ajay Melekamburath on 12/5/22.
//

#include "cc.h"

int_struct get_integrals(const scf_results& SCF, const mp2_results& MP2){
    int_struct output;
    output.oooo = slice_ints(MP2.so_ints, SCF.noo, SCF.nvo, "oooo");
    output.ovoo = slice_ints(MP2.so_ints, SCF.noo, SCF.nvo, "ovoo");
    output.oovo = slice_ints(MP2.so_ints, SCF.noo, SCF.nvo, "oovo");
    output.ooov = slice_ints(MP2.so_ints, SCF.noo, SCF.nvo, "ooov");
    output.ovov = slice_ints(MP2.so_ints, SCF.noo, SCF.nvo, "ovov");
    output.ovvo = slice_ints(MP2.so_ints, SCF.noo, SCF.nvo, "ovvo");
    output.oovv = slice_ints(MP2.so_ints, SCF.noo, SCF.nvo, "oovv");
    output.vovv = slice_ints(MP2.so_ints, SCF.noo, SCF.nvo, "vovv");
    output.ovvv = slice_ints(MP2.so_ints, SCF.noo, SCF.nvo, "ovvv");
    output.vvvo = slice_ints(MP2.so_ints, SCF.noo, SCF.nvo, "vvvo");
    output.vvvv = slice_ints(MP2.so_ints, SCF.noo, SCF.nvo, "vvvv");

    return output;
}


// MOEs in SO basis
DTensor make_F_spin_rhf(const Vector &eps, const int& nao) {
    auto n = nao * 2;
    DTensor F_spin(n, n);
    F_spin.fill(0.0);
    for (auto i = 0; i < n; i++) {
        F_spin(i, i) = eps(i / 2);
    }
    return F_spin;
}

DTensor make_F_spin_uhf(const Vector& eps_a, const Vector& eps_b, const int& nao){
    auto n = nao * 2;
    DTensor eps_so(n, n);
    for (auto i = 0; i < n; i++){
        if (i%2 == 0)
            eps_so(i, i) = eps_a(i/2);
        else if (i%2 == 1)
            eps_so(i, i) = eps_b(i/2);
    }
    return eps_so;
}

// Stanton Equation #12
DTensor make_D_ia(const DTensor& F_spin, const int& no, const int& nv){
    auto n = no + nv;
    DTensor D_ai(n, n);
    D_ai.fill(0.0);
    for (auto i = 0; i < no; i++){
        for (auto a = no; a < n; a++){
            D_ai(a, i) = F_spin(i, i) - F_spin(a, a);
        }
    }
    return D_ai;
}

// Stanton Equation #13
DTensor make_D_ijab(const DTensor& F_spin, const int& no, const int& nv){
    auto n = no + nv;
    DTensor D_abij(n,n,n,n);
    D_abij.fill(0.0);
    for (auto i = 0; i < no; i++){
        for (auto j = 0; j < no; j++){
            for (auto a = no; a < n; a++){
                for (auto b = no; b < n; b++){
                    D_abij(a, b, i, j) = F_spin(i, i) + F_spin(j, j) - F_spin(a, a) - F_spin(b, b);
                }
            }
        }
    }
    return D_abij;
}

DTensor make_tau(const DTensor& Ts, const DTensor& Td){
    auto no = Ts.extent(0);
    auto nv = Ts.extent(1);
    DTensor tau(no, no, nv, nv);
    for (auto i = 0; i < no; i++){
        for (auto j = 0; j < no; j++){
            for (auto a = 0; a < nv; a++){
                for (auto b = 0; b < nv; b++){
                    tau(i, j, a, b) = Td(i, j, a, b) + Ts(i, a) * Ts(j, b) - Ts(i, b) * Ts(j, a); // Stanton Equation #9
                }
            }
        }
    }
    return tau;
}
DTensor make_tau_bar(const DTensor& Ts, const DTensor& Td){
    auto no = Ts.extent(0);
    auto nv = Ts.extent(1);
    DTensor tau_bar(no, no, nv, nv);
    for (auto i = 0; i < no; i++){
        for (auto j = 0; j < no; j++){
            for (auto a = 0; a < nv; a++){
                for (auto b = 0; b < nv; b++){
                    tau_bar(i, j, a, b) = Td(i, j, a, b) +
                                          0.5 * (Ts(i, a) * Ts(j, b) - Ts(i, b) * Ts(j, a)); // Stanton Equation #10
                }
            }
        }
    }
    return tau_bar;
}

DTensor multiply_Ts(const DTensor& Ts){
    auto no = Ts.extent(0);
    auto nv = Ts.extent(1);
    DTensor result(no, no, nv, nv);
    for (auto j = 0; j < no; j++){
        for (auto n = 0; n < no; n++){
            for (auto f = 0; f < nv; f++){
                for (auto b = 0; b < nv; b++){
                    result(j, n, f, b) = Ts(j, f) * Ts (n, b);
                }
            }
        }
    }
    return result;
}

// Intermediates Update
cc_intermediates update_intermediates(const DTensor& Ts, const DTensor& Td, const int_struct& integrals, const DTensor& F_spin){
    using btas::contract;
    cc_intermediates intermediates;
    auto no = Ts.extent(0);
    auto nv = Ts.extent(1);
    enum {a, b, i, j, m, n, e, f};

    DTensor F_ae(nv, nv); // Stanton Equation #3
    contract(-0.5, F_spin, {m, e}, Ts, {a, m}, 1.0, F_ae, {a, e});
    contract(1.0, Ts, {m, f}, integrals.ovvv, {m, a, f, e}, 1.0, F_ae, {a, e});
    contract(-0.5, make_tau_bar(Ts, Td), {m, n, a, f}, integrals.oovv, {m, n, e, f}, 1.0, F_ae, {a, e});
    for (auto a = 0; a < nv; a++){
        for (auto e = 0; e < nv; e++){
            F_ae(a, e) += (1 - (a == e)) * F_spin(a, e); // First term of Equation #3
        }
    }
    intermediates.F_ae = F_ae; // Storing F_ae

    DTensor F_mi(nv, nv); // Stanton Equation #4
    contract(0.5, Ts, {e, i}, F_spin, {m, e}, 1.0, F_mi, {m ,i});
    contract(1.0, Ts, {n, e}, integrals.ooov, {m, n, i, e}, 1.0, F_mi, {m, i});
    contract(0.5, make_tau_bar(Ts, Td), {i, n, e, f}, integrals.oovv, {m, n, e, f}, 1.0, F_mi, {m, i});
    for (auto m = 0; m < no; m++){
        for (auto i = 0; i < no; i++){
            F_mi(m, i) += (1 - (m == i)) * F_spin(m, i); // First term of Equation #4
        }
    }
    intermediates.F_mi = F_mi;

    DTensor F_me (no, nv); // Stanton Equation #5
    contract(1.0, Ts, {n, f}, integrals.oovv, {m, n, e, f}, 1.0, F_me, {m, e});
    for (auto m = 0; m < no; m++){
        for (auto e = 0; e < nv; e++){
            F_me(m, e) += F_spin(m, e); // First term of Equation #5
        }
    }
    intermediates.F_me = F_me;

    DTensor W_mnij(no, no, no, no); // Stanton Equation #6
    contract(1.0, Ts, {j, e}, integrals.ooov, {m, n, i, e}, 1.0, W_mnij, {m, n, i, j});
    contract(-1.0, Ts, {i, e}, integrals.ooov, {m, n, j, e}, 1.0, W_mnij, {m, n, i, j});
    contract(0.25, make_tau(Ts, Td), {i, j, e, f}, integrals.oovv, {m, n, e, f}, 1.0, W_mnij, {m, n, i, j});
    W_mnij += integrals.oooo; // First term of Equation #6

    intermediates.W_mnij = W_mnij;

    DTensor W_abef(nv, nv, nv, nv); // Stanton Equation #7
    contract(-1.0, Ts, {m, b}, integrals.vovv, {a, m, e, f}, 1.0, W_abef, {a, b, e, f});
    contract(1.0, Ts, {m, a}, integrals.vovv, {b, m, e, f}, 1.0, W_abef, {a, b, e, f});
    contract(0.25, make_tau(Ts, Td), {m, n, a, b}, integrals.oovv, {m, n, e, f}, 1.0, W_abef, {a, b, e, f});
    W_abef += integrals.vvvv; // First term of Equation #7

    intermediates.W_abef = W_abef;

    DTensor W_mbej(no, nv, nv, no); // Stanton Equation #8
    contract(1.0, Ts, {j, f}, integrals.ovvv, {m, b, e, f}, 1.0, W_mbej, {m, b, e, j});
    contract(-1.0, Ts, {n, b}, integrals.oovo, {m, n, e, j}, 1.0, W_mbej, {m, b, e, j});
    contract(-0.5, Td, {j, n, f, b}, integrals.oovv, {m, n, e, f}, 1.0, W_mbej, {m, b, e, j}); // Split last term into two
    contract(-1.0, multiply_Ts(Ts), {j, n, f, b}, integrals.oovv, {m, n, e, f}, 1.0, W_mbej, {m, b, e, j});
    W_mbej += integrals.ovvo; // First term of Equation #8

    intermediates.W_mbej = W_mbej;

    return intermediates;
}

DTensor make_T1(const DTensor& Ts, const DTensor& Td, const int_struct& integrals,
                const cc_intermediates& intermediates, const DTensor& D_ia, const DTensor& F_spin){
    using btas::contract;
    auto no = Ts.extent(0);
    auto nv = Ts.extent(1);
    enum{a, i, m, n, e, f};
    DTensor tempT1(no, nv);
    // Stanton Equation #1
    contract(1.0, Ts, {i, e}, intermediates.F_ae, {a, e}, 1.0, tempT1, {i, a});
    contract(-1.0, Ts, {m, a}, intermediates.F_mi, {m, i}, 1.0, tempT1, {i, a});
    contract(1.0, Td, {i, m, a, e}, intermediates.F_me, {m, e}, 1.0, tempT1, {i, a});
    contract(-1.0, Ts, {n, f}, integrals.ovov, {n, a, i, f}, 1.0, tempT1, {i, a});
    contract(-0.5, Td, {i, m, e, f}, integrals.ovvv, {m, a, e, f}, 1.0, tempT1, {i, a});
    contract(-0.5, Td, {m, n, a, e}, integrals.oovo, {n, m, e, i}, 1.0, tempT1, {i, a});
    // For first part of equation #1
    DTensor moes_ia(no, nv);
    for (auto i = 0; i < no; i++){
        for (auto a = 0; a < nv; a++){
            moes_ia(i, a) = F_spin(i, a);
        }
    }
    tempT1 += moes_ia;

    // Since Equation #1 is for t_ia * D_ia
    DTensor newT1(no, nv);
    for (auto i = 0; i < no; i++) {
        for (auto a = 0; a < nv; a++) {
            newT1(i, a) = tempT1(i,a)/D_ia(i, a);
        }
    }
    tempT1(0, 0); // Free up
    return newT1;
}

DTensor make_T2(const DTensor& Ts, const DTensor& Td, const int_struct& integrals,
                const cc_intermediates& intermediates, const DTensor& D_ijab, const DTensor& F_spin){
    using btas::contract;
    auto no = Ts.extent(0);
    auto nv = Ts.extent(1);

    DTensor tempT2(no, no, nv, nv);
    enum{a, b, i, j, m, n, e, f};
    // Stanton Equation #2
    tempT2 += integrals.oovv; // <ij||ab>

    // P_(ab) t_ijae F_be
    contract(1.0, Td, {i, j, a, e}, intermediates.F_ae, {b, e}, 1.0, tempT2, {i, j, a, b});
    contract(-1.0, Td, {i, j, b, e}, intermediates.F_ae, {a, e}, 1.0, tempT2, {i, j, a, b});

    // P_(ab) t_ijae 0.5 * t_mb F_me
    // Do the inner contraction first
    DTensor d1(nv, nv); // Temp tensor
    contract(1.0, Ts, {m, b}, intermediates.F_me, {m, e}, 1.0, d1, {b, e});
    // Contraction for rest of the equation
    contract(-0.5, Td, {i, j, a, e}, d1, {b, e}, 1.0, tempT2, {i, j, a, b});
    // Do the same with switching a & b
    d1 -= d1; // Clearing d1
    contract(1.0, Ts, {m, a}, intermediates.F_me, {m, e}, 1.0, d1, {a, e});
    contract(-0.5, Td, {i, j, b, e}, d1, {a, e}, 1.0, tempT2, {i, j, a, b});
    d1(0, 0); // Free up

    // -P_(ij) t_imab F_be
    contract(-1.0, Td, {i, m, a, b}, intermediates.F_mi, {m, j}, 1.0, tempT2, {i, j, a, b});
    contract(1.0, Td, {j, m, a, b}, intermediates.F_mi, {m, i}, 1.0, tempT2, {i, j, a, b});

    // -P_(ij) t_imab * t_je F_me
    DTensor d2(no, no);
    contract(1.0, Ts, {j, e}, intermediates.F_me, {m, e}, 1.0, d2, {j, m});
    contract(-0.5, Td, {i, m, a, b}, d2, {j, m}, 1.0, tempT2, {i, j, a, b});
    d2 -= d2;
    contract(1.0, Ts, {i, e}, intermediates.F_me, {m, e}, 1.0, d2, {i, m});
    contract(-0.5, Td, {j, m, a, b}, d2, {i, m}, 1.0, tempT2, {i, j, a, b});

    contract(0.5, make_tau(Ts, Td), {m, n, a, b}, intermediates.W_mnij, {m, n, i, j}, 1.0, tempT2, {i, j, a, b}); //4th term
    contract(0.5, make_tau(Ts, Td), {i, j, e, f}, intermediates.W_abef, {a, b, e, f}, 1.0, tempT2, {i, j, a, b}); //5th term

    //6th term
    // First Part
    contract(1.0, Td, {i, m, a, e}, intermediates.W_mbej, {m, b, e, j}, 1.0, tempT2, {i, j, a, b});
    contract(-1.0, Td, {j ,m ,a ,e}, intermediates.W_mbej, {m, b, e, i}, 1.0, tempT2, {i, j, a, b});
    contract(-1.0, Td, {i, m, b, e}, intermediates.W_mbej, {m, a, e, j}, 1.0, tempT2, {i, j, a, b});
    contract(1.0, Td, {j, m, b, e}, intermediates.W_mbej, {m, a, e, i}, 1.0, tempT2, {i, j, a, b});
    // Second Part
    contract(1.0, multiply_Ts(Ts), {i, m, e, a}, integrals.ovov, {m, b, e, j}, 1.0, tempT2, {i, j, a, b});
    contract(-1.0, multiply_Ts(Ts), {j, m, e, a}, integrals.ovov, {m, b, e, i}, 1.0, tempT2, {i, j, a, b});
    contract(-1.0, multiply_Ts(Ts), {i, m, e, b}, integrals.ovov, {m, a, e, j}, 1.0, tempT2, {i, j, a, b});
    contract(1.0, multiply_Ts(Ts), {j, m, e, b}, integrals.ovov, {m, a, e, i}, 1.0, tempT2, {i, j, a, b});

    //7th term
    contract(1.0, Ts, {i, e}, integrals.vvvo, {a, b, e, j}, 1.0, tempT2, {i, j, a, b});
    contract(-1.0, Ts, {j, e}, integrals.vvvo, {a, b, e, i}, 1.0, tempT2, {i, j, a, b});

    //8th term
    contract(-1.0, Ts, {m, a}, integrals.ovoo, {m, b, i, j}, 1.0, tempT2, {i, j, a, b});
    contract(1.0, Ts, {m, b}, integrals.ovoo, {m, a, i, j}, 1.0, tempT2, {i, j, a, b});

    // Since Equation #2 is for t_ia * D_ijab
    DTensor newT2(no, no, nv, nv);
    for (auto i = 0; i < no; i++){
        for (auto j = 0; j < no; j++){
            for (auto a = 0; a < nv; a++){
                for (auto b = 0; b < nv; b++){
                    newT2(i, j, a, b) = tempT2(i, j, a, b)/D_ijab(i, j, a, b);
                }
            }
        }
    }
    return newT2;
}

// TDC & HFS Equation #134
real_t ccsd_energy(const DTensor& ts, const DTensor& td, const DTensor& oovv, const DTensor& F_spin){
    real_t energy = 0;
    auto no = ts.extent(0);
    auto nv = ts.extent(1);
    for (auto i = 0; i < no; i++){
        for (auto a = 0; a < nv; a++){
            energy += F_spin(i, a) * ts(i, a);
            for (auto j = 0; j < no; j++){
                for (auto b = 0; b < nv; b++){
                    energy += 0.25 * oovv(i, j, a, b) * td(i, j, a, b) + 0.5 * oovv(i, j, a, b) * ts(i, a) * ts(j, b);
                }
            }
        }
    }
    return energy;
}

// CCSD Function

real_t CCSD(const DTensor& ERI, const scf_results& scf, const mp2_results& mp2, const params& config){
    real_t energy = 0.0;
    auto integrals = get_integrals(scf, mp2);
    DTensor F_spin;
    if (config.scf == "RHF")
        F_spin = make_F_spin_rhf(scf.moes, scf.nao);
    else if (config.scf == "UHF")
        F_spin = make_F_spin_uhf(scf.moes_a, scf.moes_b, scf.nao); // Make changes for UHF

    auto D_ia = make_D_ia(F_spin, scf.noo, scf.nvo);
    auto D_ijab = make_D_ijab(F_spin, scf.noo, scf.nvo);
    auto denom = make_denom(F_spin, scf.noo, scf.nvo);

    DTensor Ts(scf.noo, scf.nvo);
    Ts.fill(0.0);

    auto Td = mp2.T; // Initial T2 guess from MP2

    // Starting iterations
    real_t E_CC_last = 0.0;
    for (auto iter = 0; iter < config.maxiter; iter++){
        if (iter == 0)
            std::cout << "Iter         E_CC        Delta(E_CC)";
        auto E_CC = ccsd_energy(Ts, Td, integrals.oovv, F_spin);
        auto Del_E_CC = E_CC - E_CC_last;
        E_CC_last = E_CC;
        printf(" %02d %20.12f %20.12f\n", iter, E_CC, Del_E_CC);
        if (Del_E_CC < config.conv){
            break ;
        }
        // Write code for updating intermediates
        cc_intermediates intermediates = update_intermediates(Ts, Td, integrals, F_spin);
        // Make T1
        Ts = make_T1(Ts, Td, integrals, intermediates, D_ia, F_spin);
        Td = make_T2(Ts, Td, integrals, intermediates, D_ijab, F_spin);
        // Makw T2
    }

    return energy;
}
