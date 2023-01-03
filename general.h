//
// Created by Ajay Melekamburath on 11/24/22.
//

#ifndef P2_GENERAL_H
#define P2_GENERAL_H

// Some general functions

#include <iostream>
#include <vector>
#include "ext/json.hpp"

// Libint Gaussian integrals library
#include <libint2.hpp>
#if !LIBINT2_CONSTEXPR_STATICS
#include <libint2/statics_definition.h>
#endif

#include "btas/btas.h"

// Typedefs
using real_t = libint2::scalar_type;
typedef btas::Tensor<double> DTensor;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;


// Structs
struct params {
    std::string inputfile;
    std::string type;
    std::string ref;
    std::string basis;
    int charge;
    int multiplicity;
    int maxiter;
    real_t scf_conv;
    real_t cc_conv;
};

std::vector<libint2::Atom> read_geometry(const std::string &filename);
void print_geometry(const std::vector<libint2::Atom> &atoms);
//params read_config(const std::string& config_file);
params read_config_json(const std::string& config_file);
size_t nbasis(const std::vector<libint2::Shell>& shells);

#endif//P2_GENERAL_H

// EOF
