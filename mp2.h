//
// Created by Ajay Melekamburath on 11/24/22.
//

#ifndef P2_MP2_H
#define P2_MP2_H

#include "hf.h"
#include "integrals.h"

//Structs
struct mp2_results{
    real_t energy;
    DTensor so_ints;
    DTensor T;
};


//DTensor make_so_moes(const Vector& eps, const int& nao);
DTensor make_so_moes(const Vector& eps_a, const Vector& eps_b, const int& nao);

real_t mp2_energy(const DTensor& oovv, const DTensor& denom);
DTensor mp2_tensor(const DTensor& oovv, const DTensor& denom);

mp2_results MP2(const libint2::BasisSet& obs, const scf_results& scf, const params& config);

#endif//P2_MP2_H

// EOF
