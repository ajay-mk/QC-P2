//
// Created by Ajay Melekamburath on 5/14/23.
//

#include "integrals.h"
#include "core.h"

namespace qc::integrals {

real_t compute_enuc(const std::vector<libint2::Atom> &atoms) {
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

Tensor compute_ao_eri(const libint2::BasisSet &obs) {
  using libint2::Engine;
  using libint2::Operator;
  using libint2::Shell;

  const auto n = qc::utils::nbasis(obs.shells());
  Tensor ao_ints(n, n, n, n);
  ao_ints.fill(0.0);

  libint2::initialize();

  Engine engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 0);

  auto shell2bf = obs.shell2bf();

  const auto &buf = engine.results();

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

                  ao_ints(bf1, bf2, bf3, bf4) = value;
                  ao_ints(bf1, bf2, bf4, bf3) = value;
                  ao_ints(bf2, bf1, bf3, bf4) = value;
                  ao_ints(bf2, bf1, bf4, bf3) = value;
                  ao_ints(bf3, bf4, bf2, bf1) = value;
                  ao_ints(bf3, bf4, bf1, bf2) = value;
                  ao_ints(bf4, bf3, bf1, bf2) = value;
                  ao_ints(bf4, bf3, bf2, bf1) = value;
                }
              }
            }
          }
        }
      }
    }
  }
  return ao_ints;
}

Tensor transform_ao_to_mo(const Tensor &eri, const Matrix &coeffs) {
  using btas::contract;
  Tensor ia_jb;
  const auto n = eri.extent(0);

  Tensor C(n, n);
  for (auto a = 0; a < n; ++a) {
    for (auto b = 0; b < n; ++b) {
      C(a, b) = coeffs(a, b);
    }
  }
  // Tensor Contractions
  Tensor pq_rl(n, n, n, n), pq_kl(n, n, n, n), pj_kl(n, n, n, n),
      ij_kl(n, n, n, n);
  const enum { p, q, r, s, i, j, k, l };  // for annotations
  contract(1.0, eri, {p, q, r, s}, C, {s, l}, 1.0, pq_rl, {p, q, r, l});
  contract(1.0, pq_rl, {p, q, r, l}, C, {r, k}, 1.0, pq_kl, {p, q, k, l});
  contract(1.0, C, {q, j}, pq_kl, {p, q, k, l}, 1.0, pj_kl, {p, j, k, l});
  contract(1.0, C, {p, i}, pj_kl, {p, j, k, l}, 1.0, ij_kl, {i, j, k, l});

  // clean up memory
  pq_kl(0, 0, 0, 0), pj_kl(0, 0, 0, 0), pq_rl(0, 0, 0, 0);
  return ij_kl;
}

Tensor transform_ao_to_mo(const Tensor &eri, const Matrix &coeff_a,
                          const Matrix &coeff_b) {
  using btas::contract;
  Tensor ia_jb;
  const int n = eri.extent(0);

  Tensor Ca(n, n);
  Tensor Cb(n, n);
  for (auto a = 0; a < n; ++a) {
    for (auto b = 0; b < n; ++b) {
      Ca(a, b) = coeff_a(a, b);
      Cb(a, b) = coeff_b(a, b);
    }
  }
  // Tensor Contractions
  Tensor pq_rl(n, n, n, n), pq_kl(n, n, n, n), pj_kl(n, n, n, n),
      ij_kl(n, n, n, n);
  enum { p, q, r, s, i, j, k, l };  // for annotations
  contract(1.0, eri, {p, q, r, s}, Ca, {s, l}, 1.0, pq_rl, {p, q, r, l});
  contract(1.0, pq_rl, {p, q, r, l}, Ca, {r, k}, 1.0, pq_kl, {p, q, k, l});
  contract(1.0, Cb, {q, j}, pq_kl, {p, q, k, l}, 1.0, pj_kl, {p, j, k, l});

  contract(1.0, Cb, {p, i}, pj_kl, {p, j, k, l}, 1.0, ij_kl, {i, j, k, l});

  // don't need other three tensors anymore
  pq_kl(0, 0, 0, 0), pj_kl(0, 0, 0, 0), pq_rl(0, 0, 0, 0);
  return ij_kl;
}

Tensor transform_to_so(const Tensor &ints_aa, const Tensor &ints_bb,
                       const Tensor &ints_ab) {
  const auto n = ints_aa.extent(0) * 2;
  Tensor result(n, n, n, n);
  for (auto i = 0; i < n; ++i) {
    for (auto j = 0; j < n; ++j) {
      for (auto k = 0; k < n; ++k) {
        for (auto l = 0; l < n; ++l) {
          if (i % 2 == 0 && k % 2 == 0 && j % 2 == 0 && l % 2 == 0)
            result(i, j, k, l) = ints_aa(i / 2, k / 2, j / 2, l / 2) -
                                 ints_aa(j / 2, k / 2, i / 2, l / 2);

          else if (i % 2 == 1 && k % 2 == 1 && j % 2 == 1 && l % 2 == 1)
            result(i, j, k, l) = ints_bb(i / 2, k / 2, j / 2, l / 2) -
                                 ints_bb(j / 2, k / 2, i / 2, l / 2);

          else if (i % 2 == 0 && k % 2 == 0 && j % 2 == 1 && l % 2 == 1)
            result(i, j, k, l) = ints_ab(i / 2, k / 2, j / 2, l / 2);

          else if (i % 2 == 1 && k % 2 == 1 && j % 2 == 0 && l % 2 == 0)
            result(i, j, k, l) = ints_ab(i / 2, k / 2, j / 2, l / 2);

          else if (i % 2 == 1 && k % 2 == 0 && j % 2 == 0 && l % 2 == 1)
            result(i, j, k, l) = -ints_ab(j / 2, k / 2, i / 2, l / 2);

          else if (i % 2 == 0 && k % 2 == 1 && j % 2 == 1 && l % 2 == 0)
            result(i, j, k, l) = -ints_ab(j / 2, k / 2, i / 2, l / 2);
        }
      }
    }
  }
  return result;
}

Tensor build_so_energies(const Matrix &eps_a, const Matrix &eps_b,
                         const int &nao) {
  const auto n = nao * 2;
  Tensor eps_so(n, n);
  for (auto i = 0; i < n; ++i) {
    if (i % 2 == 0)
      eps_so(i, i) = eps_a(i / 2);
    else if (i % 2 == 1)
      eps_so(i, i) = eps_b(i / 2);
  }
  return eps_so;
}

Tensor build_energy_denom(const Tensor &F, const int &no, const int &nv) {
  Tensor result(no, no, nv, nv);  // E_{ijab}
  for (auto i = 0; i < no; ++i) {
    for (auto j = 0; j < no; ++j) {
      for (auto a = 0; a < nv; ++a) {
        for (auto b = 0; b < nv; ++b) {
          result(i, j, a, b) =
              F(i, i) + F(j, j) - F(no + a, no + a) - F(no + b, no + b);
        }
      }
    }
  }
  return result;
}

Tensor slice_ints(const Tensor &so_ints, const int &no, const int &nv,
                  const std::string &shape) {
  // For shape
  auto n1 = (shape[0] == 'o') * no + (shape[0] == 'v') * nv;
  auto n2 = (shape[1] == 'o') * no + (shape[1] == 'v') * nv;
  auto n3 = (shape[2] == 'o') * no + (shape[2] == 'v') * nv;
  auto n4 = (shape[3] == 'o') * no + (shape[3] == 'v') * nv;

  auto m1 = (shape[0] == 'v') * no;
  auto m2 = (shape[1] == 'v') * no;
  auto m3 = (shape[2] == 'v') * no;
  auto m4 = (shape[3] == 'v') * no;

  Tensor int_slice(n1, n2, n3, n4);
  for (auto p = 0; p < n1; ++p) {
    for (auto q = 0; q < n2; ++q) {
      for (auto r = 0; r < n3; ++r) {
        for (auto s = 0; s < n4; ++s) {
          int_slice(p, q, r, s) = so_ints(p + m1, q + m2, r + m3, s + m4);
        }
      }
    }
  }
  return int_slice;
}

}  // namespace qc::integrals
