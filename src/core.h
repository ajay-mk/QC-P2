//
// Created by Ajay Melekamburath on 5/13/23.
//

#ifndef QC_CORE_H
#define QC_CORE_H

#include <iostream>
#include <nlohmann/json.hpp>
#include <vector>

// Libint Gaussian integrals library
#include <libint2.hpp>

#if !LIBINT2_CONSTEXPR_STATICS
#include <libint2/statics_definition.h>
#endif

#include "btas/btas.h"

// Typedefs
using real_t = libint2::scalar_type;
typedef btas::Tensor<real_t> Tensor;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    Matrix;

// Structs
struct input {
  std::string geom_file;
  std::string type;
  std::string ref;
  std::string basis;
  int charge;
  int multiplicity;
  int maxiter;
  real_t scf_conv;
  real_t cc_conv;
};
namespace qc {

namespace core {
/// @brief reads geometry from .xyz files
/// @param filename location and name of the input file
/// @return std::vector<libint2::Atom>
std::vector<libint2::Atom> read_geometry(const std::string &filename);

/// @brief function for reading elements from .json
/// @param json json input as stream
/// @param key key for the config item
/// @param def_value default value for the config item
/// @return the value of config item
template <typename T>
T read_json_item(const nlohmann::json &json, const std::string &key,
                 const T &def_value);

/// @brief prints the molecular geometry
/// @param atoms std::vector<libint2::Atom>, output from read_geometry
void print_geometry(const std::vector<libint2::Atom> &atoms);

/// @brief read configuration from .json file
/// @param config_file location and filename of the config file
/// @return input struct
input read_config(const std::string &config_file);

}  // namespace core

namespace utils {
/// @brief counts the number of basis functions
/// @param shells std::vector<libint2::Shell> &shells
/// @return number of basis functions
size_t nbasis(const std::vector<libint2::Shell> &shells);

}  // namespace utils

}  // namespace qc

#endif  // QC_CORE_H
