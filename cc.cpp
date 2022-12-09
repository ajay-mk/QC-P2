//
// Created by Ajay Melekamburath on 12/5/22.
//

#include "cc.h"


real_t ccsd_energy(const DTensor& ts, const DTensor& td, const DTensor& ij_ab, const DTensor F_spin){
    real_t energy = 0;
    auto no = ts.extent(0);
    auto nv = ts.extent(1);
    for (auto i = 0; i < no; i++){
        for (auto a = 0; a < nv; a++){
            energy += F_spin(i, a) * ts(i, a);
            for (auto j = 0; j < no; j++){
                for (auto b = 0; b < nv; b++){
                    energy += 0.25 * ij_ab(i, j, a, b) * td(i, j, a, b) + 0.5 * ij_ab(i, j, a, b) * ts(i, a) * ts(j, b);
                }
            }
        }
    }
    return energy;
}

