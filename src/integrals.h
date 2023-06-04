//
// Created by Ajay Melekamburath on 5/14/23.
//

#ifndef QC_INTEGRALS_H
#define QC_INTEGRALS_H

#include "core.h"

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

/// @brief computes ERIs in AO basis
/// @param obs libint2::BasisSet object
/// @return btas::Tensor of ERIs
Tensor compute_ao_eri(const libint2::BasisSet &obs);

/// @brief transforms ERIs to MO basis from AO basis
/// @param eri ERI tensor
/// @param coeffs coefficient matrix
/// @return btas::Tensor of MO integrals
Tensor transform_ao_to_mo(const Tensor &eri, const Matrix &coeffs);

/// @brief transforms ERIs to MO basis from AO basis - used for UHF reference
/// @param eri ERI tensor
/// @param coeff_a coefficient matrix
/// @param coeff_b coefficient matrix
/// @return btas::Tensor of MO integrals
Tensor transform_ao_to_mo(const Tensor &eri, const Matrix &coeff_a,
                          const Matrix &coeff_b);

/// @brief transforms MO spatial integrals spin-orbital basis
/// @param ints_aa MO integrals, alpha-alpha spin case
/// @param ints_bb MO integrals, beta-beta spin case
/// @param ints_ab MO integrals, mixed spin case
/// @return btas::Tensor of integrals in spin-orbital basis
Tensor transform_to_so(const Tensor &ints_aa, const Tensor &ints_bb,
                       const Tensor &ints_ab);

/// @brief builds E_{ij}^{ab} tensor with diagonal elements of the fock matrix
/// @param F fock matrix integrals in spin-orbital basis
/// @param no number of occupied orbitals in spin-orbital basis
/// @param nv number of virtual orbitals in spin-orbital basis
/// @return btas::Tensor of diagonal elements of the fock matrix
Tensor build_fock_tensor(const Tensor &F, const int &no, const int &nv);

/// @brief slices the integrals in specified shape
/// @param so_ints integrals in spin-orbital basis
/// @param no number of occupied orbitals in spin-orbital basis
/// @param nv number of virtual orbitals in spin-orbital basis
/// @param shape the string with required shape. eg: "oovv", "vvoo" etc.
/// @return btas::Tensor of the desired shape
Tensor slice_ints(const Tensor &so_ints, const int &no, const int &nv,
                  const std::string &shape);

}  // namespace qc::integrals

#endif  // QC_INTEGRALS_H
