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

DTensor make_fock(const Vector& eps1, const Vector& eps2, int n1, int n2){
    DTensor result(n2 - n1, n2 - n1);
    for (auto i = 0; i < n2-n1; i++){
        if (n1%2 == 0)
            result(i, i) = (i%2==0)*eps1(n1/2 + i/2) + (i%2==1)*eps2(n1/2 + i/2);
        else if (n1%2 == 1)
            result(i, i) = (i%2==0)*eps1(n1/2 + i/2 + 1) + (i%2==1)*eps2(n1/2 + i/2 + 1);
    }
    return result;
}

moes make_moe_tensors(const scf_results& scf, const params& config){
    moes result;

    if (config.scf == "RHF"){
        result.F = make_fock(scf.moes, scf.moes, 0, scf.noo + scf.nvo);
        result.F_ii = make_fock(scf.moes, scf.moes, 0, scf.noo);
        result.F_aa = make_fock(scf.moes, scf.moes, scf.noo, scf.noo + scf.nvo);
        DTensor F_ia(scf.noo, scf.nvo);
        for (auto i = 0; i < scf.noo; i++){
            for (auto a = 0; a < scf.nvo; a++){
                F_ia(i, a) = result.F(i, scf.noo + a);
            }
        }
        result.F_ia = F_ia;
    }

    if (config.scf == "UHF"){
        result.F = make_fock(scf.moes_a, scf.moes_b, 0, scf.noo + scf.nvo);
        result.F_ii = make_fock(scf.moes_a, scf.moes_b, 0, scf.noo);

        DTensor F_aa(scf.nvo, scf.nvo);
        for (auto a = 0; a < scf.nvo; a++){
            F_aa (a, a) = result.F(scf.noo + a, scf.noo + a);
        }
        result.F_aa = F_aa;

        DTensor F_ia(scf.noo, scf.nvo);
        for (auto i = 0; i < scf.noo; i++){
            for (auto a = 0; a < scf.nvo; a++){
                F_ia(i, a) = result.F(i, scf.noo + a);
            }
        }
    }
    return result;
}

// Stanton Equation #12
DTensor make_D_ia(const moes& moes){
    auto no = moes.F_ia.extent(0);
    auto nv = moes.F_ia.extent(1);
    DTensor D_ia(no, nv);
    for (auto i = 0; i < no; i++){
        for (auto a = 0; a < nv; a++){
            D_ia(i, a) = moes.F_ii(i, i) - moes.F_aa(a, a);
        }
    }
    return D_ia;
}

// Stanton Equation #13
DTensor make_D_ijab(const moes& moes){
    auto no = moes.F_ia.extent(0);
    auto nv = moes.F_ia.extent(1);
    DTensor D_ijab(no,no,nv,nv);
    D_ijab.fill(0.0);
    for (auto i = 0; i < no; i++){
        for (auto j = 0; j < no; j++){
            for (auto a = 0; a <nv ; a++){
                for (auto b = 0; b < nv; b++){
                    D_ijab(i, j, a, b) = moes.F_ii(i, i) + moes.F_ii(j, j) - moes.F_aa(a, a) - moes.F_aa(b, b);
                }
            }
        }
    }
    return D_ijab;
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
cc_intermediates update_intermediates(const DTensor& Ts, const DTensor& Td,
                                      const int_struct& integrals, const DTensor& F_spin, const moes& moes){
    using btas::contract;
    cc_intermediates intermediates;
    auto no = Ts.extent(0);
    auto nv = Ts.extent(1);
    enum {a, b, i, j, m, n, e, f};

    DTensor F_ae(nv, nv); // Stanton Equation #3

//    std::cout << Ts.extent(0) << " " << Ts.extent(1) << std::endl;
//    std::cout << moes.F_ia.extent(0) << " " << moes.F_ia.extent(1) << std::endl;
//    std::cout << integrals.ovvv.extent(0) << " " << integrals.ovvv.extent(1) << std::endl;
    contract(-0.5, moes.F_ia, {m, e}, Ts, {m, a}, 1.0, F_ae, {a, e});
    contract(1.0, Ts, {m, f}, integrals.ovvv, {m, a, f, e}, 1.0, F_ae, {a, e});
    contract(-0.5, make_tau_bar(Ts, Td), {m, n, a, f}, integrals.oovv, {m, n, e, f}, 1.0, F_ae, {a, e});

    for (auto a = 0; a < nv; a++){
        for (auto e = 0; e < nv; e++){
            F_ae(a, e) += (1 - (a == e)) * moes.F_aa(a, e); // First term of Equation #3
        }
    }
    intermediates.F_ae = F_ae; // Storing F_ae
    //std::cout << "Eqn #3 Done" << std::endl;


    DTensor F_mi(no, no); // Stanton Equation #4
    contract(0.5, Ts, {i, e}, moes.F_ia, {m, e}, 1.0, F_mi, {m ,i});
    contract(1.0, Ts, {n, e}, integrals.ooov, {m, n, i, e}, 1.0, F_mi, {m, i});
    contract(0.5, make_tau_bar(Ts, Td), {i, n, e, f}, integrals.oovv, {m, n, e, f}, 1.0, F_mi, {m, i});

    for (auto m = 0; m < no; m++){
        for (auto i = 0; i < no; i++){
            F_mi(m, i) += (1 - (m == i)) * moes.F_ii(m, i); // First term of Equation #4
        }
    }
    intermediates.F_mi = F_mi;
    //std::cout << "Eqn #4 Done" << std::endl;

    DTensor F_me (no, nv); // Stanton Equation #5
    contract(1.0, Ts, {n, f}, integrals.oovv, {m, n, e, f}, 1.0, F_me, {m, e});
    for (auto m = 0; m < no; m++){
        for (auto e = 0; e < nv; e++){
            F_me(m, e) += moes.F_ia(m, e); // First term of Equation #5
        }
    }
    intermediates.F_me = F_me;
    //std::cout << "Eqn #5 Done" << std::endl;

    DTensor W_mnij(no, no, no, no); // Stanton Equation #6
    contract(1.0, Ts, {j, e}, integrals.ooov, {m, n, i, e}, 1.0, W_mnij, {m, n, i, j});
    contract(-1.0, Ts, {i, e}, integrals.ooov, {m, n, j, e}, 1.0, W_mnij, {m, n, i, j});
    contract(0.25, make_tau(Ts, Td), {i, j, e, f}, integrals.oovv, {m, n, e, f}, 1.0, W_mnij, {m, n, i, j});
    W_mnij += integrals.oooo; // First term of Equation #6

    intermediates.W_mnij = W_mnij;
    //std::cout << "Eqn #6 Done" << std::endl;

    DTensor W_abef(nv, nv, nv, nv); // Stanton Equation #7
    contract(-1.0, Ts, {m, b}, integrals.vovv, {a, m, e, f}, 1.0, W_abef, {a, b, e, f});
    contract(1.0, Ts, {m, a}, integrals.vovv, {b, m, e, f}, 1.0, W_abef, {a, b, e, f});
    contract(0.25, make_tau(Ts, Td), {m, n, a, b}, integrals.oovv, {m, n, e, f}, 1.0, W_abef, {a, b, e, f});
    W_abef += integrals.vvvv; // First term of Equation #7

    intermediates.W_abef = W_abef;
    //std::cout << "Eqn #7 Done" << std::endl;

    DTensor W_mbej(no, nv, nv, no); // Stanton Equation #8
    contract(1.0, Ts, {j, f}, integrals.ovvv, {m, b, e, f}, 1.0, W_mbej, {m, b, e, j});
    contract(-1.0, Ts, {n, b}, integrals.oovo, {m, n, e, j}, 1.0, W_mbej, {m, b, e, j});
    contract(-0.5, Td, {j, n, f, b}, integrals.oovv, {m, n, e, f}, 1.0, W_mbej, {m, b, e, j}); // Split last term into two
    contract(-1.0, multiply_Ts(Ts), {j, n, f, b}, integrals.oovv, {m, n, e, f}, 1.0, W_mbej, {m, b, e, j});
    W_mbej += integrals.ovvo; // First term of Equation #8

    intermediates.W_mbej = W_mbej;
    //std::cout << "Eqn #8 Done" << std::endl;

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
    //tempT1(0, 0); // Free up
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
    contract(1.0, multiply_Ts(Ts), {i, m, e, a}, integrals.ovvo, {m, b, e, j}, 1.0, tempT2, {i, j, a, b});
    contract(-1.0, multiply_Ts(Ts), {j, m, e, a}, integrals.ovvo, {m, b, e, i}, 1.0, tempT2, {i, j, a, b});
    contract(-1.0, multiply_Ts(Ts), {i, m, e, b}, integrals.ovvo, {m, a, e, j}, 1.0, tempT2, {i, j, a, b});
    contract(1.0, multiply_Ts(Ts), {j, m, e, b}, integrals.ovvo, {m, a, e, i}, 1.0, tempT2, {i, j, a, b});

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
real_t ccsd_energy(const DTensor& ts, const DTensor& td, const DTensor& oovv, const moes& moes){
    real_t energy = 0;
    auto no = ts.extent(0);
    auto nv = ts.extent(1);
    for (auto i = 0; i < no; i++){
        for (auto a = 0; a < nv; a++){
            energy += moes.F_ia(i, a) * ts(i, a);
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

real_t CCSD(const scf_results& scf, const mp2_results& mp2, const params& config){
    real_t energy = 0.0;
    auto integrals = get_integrals(scf, mp2);
    DTensor F_spin;
    if (config.scf == "RHF")
        F_spin = make_F_spin_rhf(scf.moes, scf.nao);
    else if (config.scf == "UHF")
        F_spin = make_F_spin_uhf(scf.moes_a, scf.moes_b, scf.nao); // Make changes for UHF

    auto moes = make_moe_tensors(scf, config);
    auto D_ia = make_D_ia(moes);
    auto D_ijab = make_D_ijab(moes);

    DTensor Ts(scf.noo, scf.nvo);
    Ts.fill(0.0);

    auto Td = mp2.T; // Initial T2 guess from MP2

    // Starting iterations
    real_t E_CC_last = 0.0;
    for (auto iter = 0; iter < config.maxiter; iter++){
        if (iter == 0){
            std::cout << std::endl << "Iter         E_CC (Eh)        Delta(E_CC)" << std::endl;}

        auto E_CC = ccsd_energy(Ts, Td, integrals.oovv, moes);
        auto Del_E_CC = E_CC - E_CC_last;

        E_CC_last = E_CC;

        printf(" %02d %20.12f %20.12f\n", iter, E_CC, Del_E_CC);
        if (abs(Del_E_CC) < config.conv){
            break ;
        }
        // Write code for updating intermediates
        //std::cout << "Updating intermediates" << std::endl;
        cc_intermediates intermediates = update_intermediates(Ts, Td, integrals, F_spin, moes);

        // Make T1 & T2

        Ts = make_T1(Ts, Td, integrals, intermediates, D_ia, F_spin);
        //std::cout << "T1 Updated" << std::endl;
        Td = make_T2(Ts, Td, integrals, intermediates, D_ijab, F_spin);
        //std::cout << "T2 Updated" << std::endl;

    }
    return energy;
}
