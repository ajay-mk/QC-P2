//
// Created by Ajay Melekamburath on 5/14/23.
//

#ifndef QC_INTEGRALS_H
#define QC_INTEGRALS_H

#include "utils.h"

// Libint Gaussian integrals library
#include <libint2/chemistry/sto3g_atomic_density.h>
#include <libint2.hpp>

namespace qc::integrals {

/// @brief computes nuclear-repulsion
/// @param atoms std::vector<libint2::Atom>
/// @return the energy
real_t compute_enuc(const std::vector<libint2::Atom> &atoms);

/// @brief computes one-body integrals
/// @param obs libint2::BasisSet object
/// @param obtype libint2::Operator type {overlap, kinetic, nuclear}
/// @param atoms std::vector<libint2::Atom>
/// @return matrix of one-body integrals
Matrix compute_1body_ints(const libint2::BasisSet &obs,
                          libint2::Operator obtype,
                          const std::vector<libint2::Atom> &atoms);

}  // namespace qc::integrals

#endif  // QC_INTEGRALS_H
