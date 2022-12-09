//
// Created by Ajay Melekamburath on 12/5/22.
//

#ifndef P2_CC_H
#define P2_CC_H

#include "general.h"
#include "integrals.h"
#include "hf.h"
#include "mp2.h"


real_t ccsd_energy(const DTensor& T1, const DTensor&td, const DTensor& ij_ab, const DTensor F_spin);


#endif//P2_CC_H
