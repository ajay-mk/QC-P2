//
// Created by Ajay Melekamburath on 12/2/22.
//

#ifndef P2_INTEGRALS_H
#define P2_INTEGRALS_H

#include "general.h"

// Structs
struct int_struct{
    DTensor oooo;
    DTensor ovoo;
    DTensor oovo;
    DTensor ooov;
    DTensor ovov;
    DTensor ovvo;
    DTensor oovv;
    DTensor vovv;
    DTensor ovvv;
    DTensor vvvo;
    DTensor vvvv;
};




DTensor eri_ao_tensor(const libint2::BasisSet& obs);

DTensor transform_ao_mo(const DTensor& pq_rs, const Matrix& Coeff1, const Matrix& Coeff2);

DTensor make_denom(const DTensor& F_spin, int no, int nv);

DTensor transform_to_so(const DTensor& mo_ints_aa, const DTensor& mo_ints_bb, const DTensor& mo_ints_ab);

DTensor slice_ints(const DTensor& so_ints, int no, int nv, std::string int_string);

#endif//P2_INTEGRALS_H

// EOF
