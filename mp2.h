//
// Created by Ajay Melekamburath on 11/24/22.
//

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

DTensor get_iajb(const DTensor& ij_kl, const int& nocc, const int& nuocc);
real_t mp2_energy(DTensor ia_jb, Vector eps);

mp2_results MP2(const libint2::BasisSet& obs, const scf_results& scf);
