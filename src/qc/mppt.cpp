//
// Created by Ajay Melekamburath on 5/25/23.
//

#include "mppt.h"

namespace qc {
std::tuple<Tensor, real_t> MP2::run_mp2(const Tensor& oovv, const Tensor& D) {
  const auto no = oovv.extent(0);
  const auto nv = oovv.extent(2);
  real_t energy = 0.0;
  Tensor T(no, no, nv, nv);
  for (auto i = 0; i < no; ++i) {
    for (auto j = 0; j < no; ++j) {
      for (auto a = 0; a < nv; ++a) {
        for (auto b = 0; b < nv; ++b) {
          T(i, j, a, b) = oovv(i, j, a, b) / D(i, j, a, b);
          energy +=
              0.25 * (oovv(i, j, a, b) * oovv(i, j, a, b)) / D(i, j, a, b);
        }
      }
    }
  }

  return std::tie(T, energy);
}

}  // namespace qc
