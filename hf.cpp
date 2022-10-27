// Contains functions relevant to Hartree Fock Algorithm
//
// Created by Ajay Melekamburath on 10/26/22.
//

#include <iostream>
#include <string>
#include <vector>


#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

// Libint Gaussian integrals library
#include <libint2.hpp>
#if !LIBINT2_CONSTEXPR_STATICS
#  include <libint2/statics_definition.h>
#endif

using real_t = libint2::scalar_type;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
