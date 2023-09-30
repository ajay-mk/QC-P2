//
// Created by Ajay Melekamburath on 5/13/23.
//

#include <qc/scf.h>

namespace qc {

// SCF Class functions

Matrix SCF::density_guess(int nocc, int nao) {
  Matrix result = Matrix::Zero(nao, nao);
  for (int i = 0; i < nocc; ++i) result(i, i) = 1.0;
  return result;
}

Matrix SCF::compute_soad(const std::vector<libint2::Atom> &atoms) {
  std::cout << "Using SOAD as initial guess" << std::endl;
  // compute number of atomic orbitals
  size_t nao = 0;
  for (const auto &atom : atoms) {
    const auto Z = atom.atomic_number;
    nao += libint2::sto3g_num_ao(Z);
  }
  // compute the minimal basis density
  Matrix D = Matrix::Zero(nao, nao);
  size_t ao_offset = 0;  // first AO of this atom
  for (const auto &atom : atoms) {
    const auto Z = atom.atomic_number;
    const auto &occvec = libint2::sto3g_ao_occupation_vector(Z);
    for (const auto &occ : occvec) {
      D(ao_offset, ao_offset) = occ;
      ++ao_offset;
    }
  }
  return D * 0.5;
}

Matrix SCF::build_rhf_fock(const libint2::BasisSet &obs, const Matrix &D) {
  using libint2::Engine;
  using libint2::Operator;
  using libint2::Shell;

  const auto n = qc::utils::nbasis(obs.shells());
  Matrix G = Matrix::Zero(n, n);

  // construct the 2-electron repulsion integrals engine
  Engine engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 0);

  auto shell2bf = obs.shell2bf();

  const auto &buf = engine.results();

  // loop over permutationally-unique set of shells
  for (auto s1 = 0; s1 != obs.shells().size(); ++s1) {
    auto bf1_first = shell2bf[s1];  // first basis function in this shell
    auto n1 =
        obs.shells()[s1].size();  // number of basis functions in this shell

    for (auto s2 = 0; s2 <= s1; ++s2) {
      auto bf2_first = shell2bf[s2];
      auto n2 = obs.shells()[s2].size();

      for (auto s3 = 0; s3 <= s1; ++s3) {
        auto bf3_first = shell2bf[s3];
        auto n3 = obs.shells()[s3].size();

        const auto s4_max = (s1 == s3) ? s2 : s3;
        for (auto s4 = 0; s4 <= s4_max; ++s4) {
          auto bf4_first = shell2bf[s4];
          auto n4 = obs.shells()[s4].size();

          // compute the permutational degeneracy (i.e. # of equivalents) of the
          // given shell set
          auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
          auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
          auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
          auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

          engine.compute(obs.shells()[s1], obs.shells()[s2], obs.shells()[s3],
                         obs.shells()[s4]);
          const auto *buf_1234 = buf[0];
          if (buf_1234 == nullptr)
            continue;  // if all integrals screened out, skip to next quartet

          for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for (auto f2 = 0; f2 != n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for (auto f3 = 0; f3 != n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                  const auto bf4 = f4 + bf4_first;

                  const auto value = buf_1234[f1234];

                  const auto value_scal_by_deg = value * s1234_deg;

                  G(bf1, bf2) += D(bf3, bf4) * value_scal_by_deg;
                  G(bf3, bf4) += D(bf1, bf2) * value_scal_by_deg;
                  G(bf1, bf3) -= 0.25 * D(bf2, bf4) * value_scal_by_deg;
                  G(bf2, bf4) -= 0.25 * D(bf1, bf3) * value_scal_by_deg;
                  G(bf1, bf4) -= 0.25 * D(bf2, bf3) * value_scal_by_deg;
                  G(bf2, bf3) -= 0.25 * D(bf1, bf4) * value_scal_by_deg;
                }
              }
            }
          }
        }
      }
    }
  }

  // symmetrize the result and return
  Matrix Gt = G.transpose();
  return 0.5 * (G + Gt);
}

Matrix SCF::build_uhf_fock(const libint2::BasisSet &obs, const Matrix &D,
                           const Matrix &Ds) {
  using libint2::Engine;
  using libint2::Operator;
  using libint2::Shell;

  const auto n = qc::utils::nbasis(obs.shells());
  Matrix G = Matrix::Zero(n, n);

  // construct the 2-electron repulsion integrals engine
  Engine engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 0);

  auto shell2bf = obs.shell2bf();

  const auto &buf = engine.results();

  // loop over permutationally-unique set of shells
  for (auto s1 = 0; s1 != obs.shells().size(); ++s1) {
    auto bf1_first = shell2bf[s1];  // first basis function in this shell
    auto n1 =
        obs.shells()[s1].size();  // number of basis functions in this shell

    for (auto s2 = 0; s2 <= s1; ++s2) {
      auto bf2_first = shell2bf[s2];
      auto n2 = obs.shells()[s2].size();

      for (auto s3 = 0; s3 <= s1; ++s3) {
        auto bf3_first = shell2bf[s3];
        auto n3 = obs.shells()[s3].size();

        const auto s4_max = (s1 == s3) ? s2 : s3;
        for (auto s4 = 0; s4 <= s4_max; ++s4) {
          auto bf4_first = shell2bf[s4];
          auto n4 = obs.shells()[s4].size();

          // compute the permutational degeneracy (i.e. # of equivalents) of the
          // given shell set
          auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
          auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
          auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
          auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

          engine.compute(obs.shells()[s1], obs.shells()[s2], obs.shells()[s3],
                         obs.shells()[s4]);
          const auto *buf_1234 = buf[0];
          if (buf_1234 == nullptr)
            continue;  // if all integrals screened out, skip to next quartet

          for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for (auto f2 = 0; f2 != n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for (auto f3 = 0; f3 != n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                  const auto bf4 = f4 + bf4_first;

                  const auto value = buf_1234[f1234];

                  const auto value_scal_by_deg = value * s1234_deg;

                  G(bf1, bf2) += 0.5 * D(bf3, bf4) * value_scal_by_deg;
                  G(bf3, bf4) += 0.5 * D(bf1, bf2) * value_scal_by_deg;
                  G(bf1, bf3) -= 0.25 * Ds(bf2, bf4) * value_scal_by_deg;
                  G(bf2, bf4) -= 0.25 * Ds(bf1, bf3) * value_scal_by_deg;
                  G(bf1, bf4) -= 0.25 * Ds(bf2, bf3) * value_scal_by_deg;
                  G(bf2, bf3) -= 0.25 * Ds(bf1, bf4) * value_scal_by_deg;
                }
              }
            }
          }
        }
      }
    }
  }

  // symmetrize the result and return
  Matrix Gt = G.transpose();
  return 0.5 * (G + Gt);
}

real_t SCF::rhf_energy(const Matrix &D, const Matrix &H, const Matrix &F) {
  real_t energy = 0.0;
  for (auto i = 0; i < D.rows(); ++i)
    for (auto j = 0; j < D.rows(); ++j) energy += D(i, j) * (H(i, j) + F(i, j));
  return energy;
}

real_t SCF::uhf_energy(const Matrix &D, const Matrix &Dalpha,
                       const Matrix &Dbeta, const Matrix &H,
                       const Matrix &Falpha, const Matrix &Fbeta) {
  real_t energy = 0.0;
  for (auto i = 0; i < D.rows(); ++i)
    for (auto j = 0; j < D.rows(); ++j)
      energy += D(i, j) * H(i, j) + Dalpha(i, j) * Falpha(i, j) +
                Dbeta(i, j) * Fbeta(i, j);
  return 0.5 * energy;
}

}  // namespace qc
