//
// Created by Ajay Melekamburath on 5/25/23.
//

#ifndef QC_MPPT_H
#define QC_MPPT_H

#include <tuple>
#include "core.h"
#include "scf.h"

namespace qc {

class MP2 {
 public:
  Tensor mp2_amps;
  real_t corr_energy;

  MP2(const std::vector<libint2::Atom>& atoms, const input& config) {
    const auto scf = SCF(atoms, config);

    std::cout << "Starting MP2 calculation\n" << std::endl;
    const auto no = scf.nocc;
    const auto nv = scf.nvir;
    const auto nao = scf.nao;

    const auto ao_ints = integrals::compute_ao_eri(scf.basis);

    if (config.ref == "RHF") {
      const auto mo_ints = integrals::transform_ao_to_mo(ao_ints, scf.C);
      const auto so_ints =
          integrals::transform_to_so(mo_ints, mo_ints, mo_ints);
      const auto F = integrals::build_so_energies(scf.eps, scf.eps, nao);
      const auto D = integrals::build_energy_denom(F, no, nv);
      const auto oovv = integrals::slice_ints(so_ints, no, nv, "oovv");

      std::tie(mp2_amps, corr_energy) = run_mp2(oovv, D);
    }  // enf of RHF block
    else if (config.ref == "UHF") {
      const auto mo_ints_aa =
          integrals::transform_ao_to_mo(ao_ints, scf.Calpha, scf.Calpha);
      const auto mo_ints_bb =
          integrals::transform_ao_to_mo(ao_ints, scf.Cbeta, scf.Cbeta);
      const auto mo_ints_ab =
          integrals::transform_ao_to_mo(ao_ints, scf.Calpha, scf.Cbeta);

      const auto so_ints =
          integrals::transform_to_so(mo_ints_aa, mo_ints_bb, mo_ints_ab);

      const auto F =
          integrals::build_so_energies(scf.eps_alpha, scf.eps_beta, nao);
      const auto D = integrals::build_energy_denom(F, no, nv);
      const auto oovv = integrals::slice_ints(so_ints, no, nv, "oovv");

      std::tie(mp2_amps, corr_energy) = run_mp2(oovv, D);
    }
    std::cout << "MP2 correlation energy = " << corr_energy << std::endl;
  }

  /// @brief computes MP2 amplitudes and correlation energy
  /// @param oovv slice of spin-orbital integrals in ijab shape
  /// @param D energy denominator in spin-orbital basis
  /// returns std::tuple of MP2 amplitudes and correlation energy
  std::tuple<Tensor, real_t> run_mp2(const Tensor& oovv, const Tensor& D);

};  // Class MP2

}  // namespace qc

#endif  // QC_MPPT_H
