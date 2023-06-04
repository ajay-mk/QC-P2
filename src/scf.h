//
// Created by Ajay Melekamburath on 5/13/23.
//

#ifndef QC_SCF_H
#define QC_SCF_H

#include <iostream>
#include <istream>
#include <string>
#include <vector>

#include <Eigen/Eigenvalues>

#include "core.h"
#include "integrals.h"

namespace qc {

class SCF {
 private:
  Matrix S, T, V, H;

 protected:
  int nelectron, nalpha, nbeta;
  int nao, nocc, nvir;  // these are in so basis
  libint2::BasisSet basis;

 public:
  real_t scf_energy;
  Matrix D, Dalpha, Dbeta;
  Matrix F, Falpha, Fbeta;
  Matrix C, Calpha, Cbeta;
  Matrix eps, eps_alpha, eps_beta;
  real_t nuclear_repulsion;

  /// @brief computes scf energy for RHF or UHF reference
  /// @param atoms std::vector<libint2::Atom>
  /// @param config config json object
  SCF(const std::vector<libint2::Atom> &atoms, const input &config) {
    core::print_geometry(atoms);
    // Counting the number of electrons
    nelectron = 0;
    for (auto &atom : atoms) nelectron += atom.atomic_number;
    std::cout << "\nNumber of electrons = " << nelectron << std::endl;
    auto ndocc = nelectron / 2;
    nuclear_repulsion = integrals::compute_enuc(atoms);

    using libint2::BasisSet;
    const auto verbose = config.verbose;  // additional printing

    libint2::BasisSet basis(config.basis, atoms);
    nao = qc::utils::nbasis(basis.shells());
    std::cout << "Number of basis functions = " << nao << std::endl;

    // Occupied and Virtual Orbitals
    nocc = 2 * (nelectron / 2);
    nvir = (2 * nao) - nocc;
    std::cout << "Number of occupied orbitals: " << nocc << std::endl
              << "Number of virtual orbitals: " << nvir << std::endl;

    std::cout << "\nStarting SCF calculation" << std::endl;
    scf_energy = 0.0;  // initialize scf energy;

    // Initializing Libint
    libint2::initialize();

    S = integrals::compute_1body_ints(basis, libint2::Operator::overlap, atoms);
    T = integrals::compute_1body_ints(basis, libint2::Operator::kinetic, atoms);
    V = integrals::compute_1body_ints(basis, libint2::Operator::nuclear, atoms);
    if (verbose) {
      std::cout << "Overlap matrix: \n" << S << std::endl;
      std::cout << "Kinetic energy matrix: \n" << T << std::endl;
      std::cout << "Potential energy matrix: \n" << V << std::endl;
    }

    H = T + V;  // Core Hamiltonian = T + V
    if (verbose) {
      std::cout << "Core Hamiltonian: \n" << H << std::endl;
    }

    // T and V no longer needed, free up the memory
    T.resize(0, 0);
    V.resize(0, 0);

    // from here, split the rhf and uhf stuff
    if (config.ref == "RHF") {
      std::cout << "Reference: " << config.ref << std::endl;

      D = (config.basis == "STO-3G" || config.basis == "sto-3g")
              ? compute_soad(atoms)
              : density_guess(nocc, nao);
      // SCF Loop
      real_t rmsd, ediff, ehf;

      for (auto iter = 0; iter < config.maxiter;) {
        ++iter;
        auto ehf_last = ehf;
        auto D_last = D;

        F = H + build_rhf_fock(basis.shells(), D);
        if (verbose && iter == 1) {
          std::cout << "Initial Fock Matrix: \n" << F << std::endl;
        }

        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver(F, S);
        eps = solver.eigenvalues();
        C = solver.eigenvectors();

        auto C_occ = C.leftCols(ndocc);
        D = C_occ * C_occ.transpose();

        ehf = rhf_energy(D, H, F);

        ediff = ehf - ehf_last;
        rmsd = (D - D_last).norm();

        if (iter == 1)
          std::cout << "\nIter\tE(elec)\tE(tot)\tDelta(E)\tRMS(D)\n";
        printf(" %02d %20.12f %20.12f %20.12f %20.12f\n", iter, ehf,
               ehf + nuclear_repulsion, ediff, rmsd);

        if (fabs(ediff) < config.scf_conv && fabs(rmsd) < config.scf_conv) {
          break;
        } else
          continue;
        if (iter >= config.maxiter) std::cout << "Iter > MaxIter" << std::endl;
      }
      scf_energy = ehf + nuclear_repulsion;
    }  // RHF block

    if (config.ref == "UHF") {
      std::cout << "Reference: " << config.ref << std::endl;
      nbeta = (nelectron - config.multiplicity + 1) / 2;
      nalpha = nbeta + config.multiplicity - 1;

      std::cout << "Number of alpha electrons: " << nalpha << std::endl
                << "Number of beta electrons: " << nbeta << std::endl;

      Dalpha = density_guess(nalpha, nao);
      Dbeta = density_guess(nbeta, nao);
      D = Dalpha + Dbeta;  // Total Density Matrix
      // SCF Loop
      real_t rmsd, ediff, euhf;

      for (auto iter = 0; iter < config.maxiter;) {
        ++iter;
        auto euhf_last = euhf;
        auto D_last = D;

        // New Fock Matrices
        Falpha = H;
        Falpha += build_uhf_fock(basis.shells(), D, Dalpha);
        Fbeta = H;
        Fbeta += build_uhf_fock(basis.shells(), D, Dbeta);

        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver1(Falpha, S);
        eps_alpha = solver1.eigenvalues();
        Calpha = solver1.eigenvectors();

        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver2(Fbeta, S);
        eps_beta = solver2.eigenvalues();
        Cbeta = solver2.eigenvectors();

        // Density Matrices
        auto Ca_occ = Calpha.leftCols(nalpha);
        Dalpha = Ca_occ * Ca_occ.transpose();

        auto Cb_occ = Cbeta.leftCols(nbeta);
        Dbeta = Cb_occ * Cb_occ.transpose();

        D = Dalpha + Dbeta;

        // UHF Energy
        euhf = uhf_energy(D, Dalpha, Dbeta, H, Falpha, Fbeta);

        // compute difference with last iteration
        ediff = euhf - euhf_last;
        rmsd = (D - D_last).norm();

        if (iter == 1)
          std::cout << "\n\n Iter        E(elec)              E(tot)           "
                       "    Delta(E)             RMS(D)\n";
        printf(" %02d %20.12f %20.12f %20.12f %20.12f\n", iter, euhf,
               euhf + nuclear_repulsion, ediff, rmsd);
        if (fabs(ediff) < config.scf_conv && fabs(rmsd) < config.scf_conv) {
          break;
        } else
          continue;
        if (iter >= config.maxiter) std::cout << "Iter > MaxIter" << std::endl;
      }
      scf_energy = euhf + nuclear_repulsion;
    }  // UHF block

    std::cout << "SCF Energy: " << scf_energy << std::endl;
  };

  // Function declarations
  /// @brief Computes Superposition-Of-Atomic-Densities guess for the molecular
  /// density matrix in minimal basis; occupies subshells by smearing electrons
  /// evenly over the orbitals
  /// @param atoms std::vector<libint2::Atom>
  /// @return matrix of initial density guess
  Matrix compute_soad(const std::vector<libint2::Atom> &atoms);

  /// @brief guess for initial density
  /// adds 1 as diagonal elements for all occupied electrons
  /// @param nocc number of occupied orbitals
  /// @param nao number of atomic orbitals
  static static Matrix density_guess(int nocc, int nao);

  /// @brief computes fock matrix for RHF reference
  /// @param obs libint2::BasisSet object
  /// @param D density matrix
  Matrix build_rhf_fock(const libint2::BasisSet &obs, const Matrix &D);

  /// @brief computes fock matrix for UHF reference
  /// @param obs libint2::BasisSet object
  /// @param D density matrix
  Matrix build_uhf_fock(const libint2::BasisSet &obs, const Matrix &D,
                        const Matrix &Ds);

  /// @brief computes RHF energy
  /// @param D density matrix
  /// @param H Hamiltonian matrix
  /// @param F fock matrix
  /// @return returns energy
  real_t rhf_energy(const Matrix &D, const Matrix &H, const Matrix &F);

  /// @brief computes UHF energy
  /// @param D initial density matrix
  /// @param Dalpha density matrix - alpha spin
  /// @param Dbeta density matrix - beta spin
  /// @param H Hamiltonian matrix
  /// @param Falpha fock matrix - alpha spin
  /// @param Fbeta fock matrix - beta spin
  /// @return returns energy
  real_t uhf_energy(const Matrix &D, const Matrix &Dalpha, const Matrix &Dbeta,
                    const Matrix &H, const Matrix &Falpha, const Matrix &Fbeta);

};  // class SCF

}  // namespace qc

#endif  // QC_SCF_H
