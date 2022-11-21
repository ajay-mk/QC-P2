//
// Created by Ajay Melekamburath on 11/14/22.
//
// Contains functions relevant to MP2 Calculations

// Libint Gaussian integrals library
#include <libint2.hpp>
#include <libint2/chemistry/sto3g_atomic_density.h>
#if !LIBINT2_CONSTEXPR_STATICS
#  include <libint2/statics_definition.h>
#endif
#include "btas/btas.h"

//TypeDefs
using real_t = libint2::scalar_type;
typedef btas::Tensor<double> Tensor;

// Functions
Tensor eri_ao_tensor(const libint2::BasisSet& obs);
size_t nbasis(const std::vector<libint2::Shell>& shells);

// Function Definitions

Tensor eri_ao_tensor(const libint2::BasisSet& obs) {
    using libint2::Shell;
    using libint2::Engine;
    using libint2::Operator;

    const auto n = nbasis(obs.shells());
    Tensor ao_ints(n, n, n, n);

    libint2::initialize();

    Engine engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 0);

    auto shell2bf = obs.shell2bf();

    const auto &buf = engine.results();

    for (auto s1 = 0; s1 != obs.shells().size(); ++s1) {

        auto bf1_first = shell2bf[s1];    // first basis function in this shell
        auto n1 = obs.shells()[s1].size();// number of basis functions in this shell

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

                    // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
                    auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
                    auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
                    auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
                    auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

                    engine.compute(obs.shells()[s1], obs.shells()[s2], obs.shells()[s3], obs.shells()[s4]);
                    const auto *buf_1234 = buf[0];
                    if (buf_1234 == nullptr)
                        continue;// if all integrals screened out, skip to next quartet

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

// EOF