//
// Created by Ajay Melekamburath on 5/14/23.
//

#include "integrals.h"

namespace qc::integrals {

double compute_enuc(const std::vector<libint2::Atom> &atoms) {
  auto num = 0.0;
  for (auto i = 0; i < atoms.size(); ++i)
    for (auto j = i + 1; j < atoms.size(); ++j) {
      auto xij = atoms[i].x - atoms[j].x;
      auto yij = atoms[i].y - atoms[j].y;
      auto zij = atoms[i].z - atoms[j].z;
      auto r2 = xij * xij + yij * yij + zij * zij;
      auto r = sqrt(r2);
      num += atoms[i].atomic_number * atoms[j].atomic_number / r;
    }
  return num;
}

// Integral engines and fock builder are from:
// https://github.com/evaleev/libint/blob/master/tests/hartree-fock/hartree-fock.cc
Matrix compute_1body_ints(const libint2::BasisSet &obs,
                          libint2::Operator obtype,
                          const std::vector<libint2::Atom> &atoms) {
  using libint2::Engine;
  using libint2::Operator;
  using libint2::Shell;

  const auto n = qc::utils::nbasis(obs.shells());
  Matrix result(n, n);

  // construct the overlap integrals engine
  Engine engine(obtype, obs.max_nprim(), obs.max_l(), 0);
  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical
  // charges
  if (obtype == Operator::nuclear) {
    std::vector<std::pair<real_t, std::array<real_t, 3>>> q;
    for (const auto &atom : atoms) {
      q.push_back({static_cast<real_t>(atom.atomic_number),
                   {{atom.x, atom.y, atom.z}}});
    }
    engine.set_params(q);
  }

  auto shell2bf = obs.shell2bf();

  // buf[0] points to the target shell set after every call  to engine.compute()
  const auto &buf = engine.results();

  // loop over unique shell pairs, {s1,s2} such that s1 >= s2
  // this is due to the permutational symmetry of the real integrals over
  // Hermitian operators: (1|2) = (2|1)
  for (auto s1 = 0; s1 != obs.shells().size(); ++s1) {
    auto bf1 = shell2bf[s1];  // first basis function in this shell
    auto n1 = obs.shells()[s1].size();

    for (auto s2 = 0; s2 <= s1; ++s2) {
      auto bf2 = shell2bf[s2];
      auto n2 = obs.shells()[s2].size();

      // compute shell pair
      engine.compute(obs.shells()[s1], obs.shells()[s2]);

      // "map" buffer to a const Eigen Matrix, and copy it to the corresponding
      // blocks of the result
      Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2)  // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1}
                     // block, note the transpose!
        result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
    }
  }
  return result;
}

}  // namespace qc::integrals
