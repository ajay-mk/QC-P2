//
// Created by Ajay Melekamburath on 12/5/22.
//

#ifndef P2_CC_H
#define P2_CC_H

#include "general.h"
#include "integrals.h"
#include "hf.h"
#include "mp2.h"

// Structs
struct cc_intermediates{
    DTensor F_ae, F_mi, F_me, W_mnij, W_mbej, W_abef;
};
struct moes{
    DTensor F, F_ii, F_ia, F_aa;
};


DTensor make_F_spin_rhf(const Vector &eps, const int& nao);
DTensor make_F_spin_uhf(const Vector &eps_a, const Vector& eps_b, const int& nao);

moes make_moe_tensors(const scf_results& scf, const params& config);

DTensor make_fock(const Vector& eps1, const Vector& eps2, int n1, int n2);

DTensor make_D_ia(const moes& moes);
DTensor make_D_ijab(const DTensor& F_spin, const int& no, const int& nv);

DTensor make_tau(const DTensor& Ts, const DTensor& Td);
DTensor make_tau_bar(const DTensor& Ts, const DTensor& Td);
DTensor multiply_Ts(const DTensor& Ts);

DTensor make_T1(const DTensor& Ts, const DTensor& Td, const int_struct& integrals,
                const cc_intermediates& intermediates, const DTensor& D_ia, const DTensor& F_spin);
DTensor make_T2(const DTensor& Ts, const DTensor& Td, const int_struct& integrals,
                const cc_intermediates& intermediates, const DTensor& D_ijab, const DTensor& F_spin);

cc_intermediates update_intermediates(const DTensor& Ts, const DTensor& Td, const int_struct& integrals, const DTensor& F_spin, const moes& moes);

real_t ccsd_energy(const DTensor& ts, const DTensor& td, const DTensor& oovv, const moes& moes);

real_t CCSD(const scf_results& SCF, const mp2_results& MP2, const params& config);

#endif//P2_CC_H
