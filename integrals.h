//
// Created by Ajay Melekamburath on 12/2/22.
//

#ifndef P2_INTEGRALS_H
#define P2_INTEGRALS_H

#include "general.h"


DTensor eri_ao_tensor(const libint2::BasisSet& obs);

DTensor transform_ao_mo(const DTensor& pq_rs, const Matrix& C);
DTensor transform_ao_mo_uhf(const DTensor& pq_rs, const Matrix& Coeff1, const Matrix& Coeff2);

DTensor transform_to_so(const DTensor& mo_ints);
DTensor transform_to_so_uhf(const DTensor& mo_ints_aa, const DTensor& mo_ints_bb, const DTensor& mo_ints_ab);

#endif//P2_INTEGRALS_H
