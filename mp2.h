//
// Created by Ajay Melekamburath on 11/24/22.
//

#ifndef P2_MP2_H
#define P2_MP2_H

#include "hf.h"

//TypeDefs
typedef btas::Tensor<double> DTensor;

//Structs
struct mp2_results{
    real_t energy;
    DTensor T;
};

DTensor transform_ao_mo(const DTensor& pq_rs, const Matrix& C);
DTensor eri_ao_tensor(const libint2::BasisSet& obs);
DTensor transform_ao_mo(const DTensor& ao_ints, const Matrix& C, const int& nocc, const int& nuocc);

Vector make_so_moes(const Vector& eps);
DTensor transform_to_so(const DTensor& mo_ints);
DTensor get_ijab(const DTensor& so_ints, const int& noo, const int& nvo);
real_t mp2_energy(DTensor ij_ab, Vector eps_so);

//DTensor get_iajb(const DTensor& ij_kl, const int& nocc, const int& nuocc);
//real_t mp2_energy(DTensor ia_jb, Vector eps);

mp2_results MP2(const libint2::BasisSet& obs, const scf_results& scf, const params& config);


#endif//P2_MP2_H

// EOF
