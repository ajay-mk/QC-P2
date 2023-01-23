//
// Created by Ajay Melekamburath on 11/24/22.
//
#ifndef P2_HF_H
#define P2_HF_H

#include <iostream>
#include <istream>
#include <string>
#include <vector>

#include <Eigen/Eigenvalues>

// Libint Gaussian integrals library
#include <libint2.hpp>
#include <libint2/chemistry/sto3g_atomic_density.h>
#if !LIBINT2_CONSTEXPR_STATICS
#include <libint2/statics_definition.h>
#endif

#include "general.h"

//TypeDefs
using real_t = libint2::scalar_type;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, 1> Vector;

struct scf_results {
    real_t energy;
    int nalpha, nbeta, no, nv;
    int nao;
    Matrix F, Fa, Fb, C, Ca, Cb, D, Da, Db;
    Vector moes, moes_a, moes_b;
};

Matrix compute_soad(const std::vector<libint2::Atom> &atoms);
double compute_enuc(const std::vector<libint2::Atom> &atoms);

Matrix compute_1body_ints(const libint2::BasisSet &obs, libint2::Operator t, const std::vector<libint2::Atom> &atoms = std::vector<libint2::Atom>());

Matrix density_guess(int nocc, int nao);
Matrix build_fock(const libint2::BasisSet &obs, const Matrix &D);
Matrix build_uhf_fock(const libint2::BasisSet &obs, const Matrix &D, const Matrix &Ds);
real_t rhf_energy(const Matrix &D, const Matrix &H, const Matrix &F);
real_t uhf_energy(const Matrix &D, const Matrix &Dalpha, const Matrix &Dbeta, const Matrix &H, const Matrix &Falpha, const Matrix &Fbeta);


scf_results RHF(const std::vector<libint2::Atom> &atoms, const libint2::BasisSet &obs, real_t nelectron, params config);
scf_results UHF(const std::vector<libint2::Atom> &atoms, const libint2::BasisSet &obs, real_t nelectron, params config);

#endif//P2_HF_H

// EOF
