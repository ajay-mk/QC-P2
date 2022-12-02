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
    DTensor T;
};


Vector make_so_moes(const Vector& eps, const int& nao);
Vector make_so_moes_uhf(const Vector& eps_a, const Vector& eps_b, const int& nao);

DTensor get_ijab(const DTensor& so_ints, const int& noo, const int& nvo);
real_t mp2_energy(DTensor ij_ab, Vector eps_so);

//DTensor get_iajb(const DTensor& ij_kl, const int& nocc, const int& nuocc);
//real_t mp2_energy(DTensor ia_jb, Vector eps);

mp2_results MP2(const libint2::BasisSet& obs, const scf_results& scf, const params& config);


#endif//P2_MP2_H

// EOF
