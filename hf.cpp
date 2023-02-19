// Contains functions relevant to Hartree-Fock Algorithm
//
// Created by Ajay Melekamburath on 10/26/22.
//

#include <iostream>
#include <istream>
#include <string>
#include <vector>

#include <Eigen/Eigenvalues>

// Include Headers
#include "hf.h"


// Function Definitions


// Computing Nuclear Repulsion Energy
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

// Integral engines are from the example: https://github.com/evaleev/libint/blob/master/tests/hartree-fock/hartree-fock.cc

Matrix
compute_1body_ints(const libint2::BasisSet &obs, libint2::Operator obtype, const std::vector<libint2::Atom> &atoms) {
    using libint2::Engine;
    using libint2::Operator;
    using libint2::Shell;

    const auto n = nbasis(obs.shells());
    Matrix result(n, n);

    // construct the overlap integrals engine
    Engine engine(obtype, obs.max_nprim(), obs.max_l(), 0);
    // nuclear attraction ints engine needs to know where the charges sit ...
    // the nuclei are charges in this case; in QM/MM there will also be classical charges
    if (obtype == Operator::nuclear) {
        std::vector<std::pair<real_t, std::array<real_t, 3>>> q;
        for (const auto &atom: atoms) {
            q.push_back({static_cast<real_t>(atom.atomic_number), {{atom.x, atom.y, atom.z}}});
        }
        engine.set_params(q);
    }

    auto shell2bf = obs.shell2bf();

    // buf[0] points to the target shell set after every call  to engine.compute()
    const auto &buf = engine.results();

    // loop over unique shell pairs, {s1,s2} such that s1 >= s2
    // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
    for (auto s1 = 0; s1 != obs.shells().size(); ++s1) {

        auto bf1 = shell2bf[s1];// first basis function in this shell
        auto n1 = obs.shells()[s1].size();

        for (auto s2 = 0; s2 <= s1; ++s2) {

            auto bf2 = shell2bf[s2];
            auto n2 = obs.shells()[s2].size();

            // compute shell pair
            engine.compute(obs.shells()[s1], obs.shells()[s2]);

            // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
            Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
            result.block(bf1, bf2, n1, n2) = buf_mat;
            if (s1 != s2)// if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
                result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
        }
    }
    return result;
}


//Computes Superposition-Of-Atomic-Densities guess for the molecular density matrix
//in minimal basis; occupies subshells by smearing electrons evenly over the orbitals
Matrix compute_soad(const std::vector<libint2::Atom> &atoms) {
    // compute number of atomic orbitals
    size_t nao = 0;
    for (const auto &atom: atoms) {
        const auto Z = atom.atomic_number;
        nao += libint2::sto3g_num_ao(Z);
    }

    // compute the minimal basis density
    Matrix D = Matrix::Zero(nao, nao);
    size_t ao_offset = 0;// first AO of this atom
    for (const auto &atom: atoms) {
        const auto Z = atom.atomic_number;
        const auto &occvec = libint2::sto3g_ao_occupation_vector(Z);
        for (const auto &occ: occvec) {
            D(ao_offset, ao_offset) = occ;
            ++ao_offset;
        }
    }

    return D * 0.5;// we use densities normalized to # of electrons/2
}
// SAD guess only works for STO-3G now, should fix this

// Guessing Initial Density - Adds 1 as diagonal elements for all occupied electrons
Matrix density_guess(int nocc, int nao) {
    Matrix guess = Matrix::Zero(nao, nao);
    for (int i = 0; i < nocc; ++i)
        guess(i, i) = 1.0;
    return guess;
}

// Fock Builder
Matrix build_fock(const libint2::BasisSet &obs, const Matrix &D) {

    using libint2::Engine;
    using libint2::Operator;
    using libint2::Shell;

    const auto n = nbasis(obs.shells());
    Matrix G = Matrix::Zero(n, n);

    // construct the 2-electron repulsion integrals engine
    Engine engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 0);

    auto shell2bf = obs.shell2bf();

    const auto &buf = engine.results();

    // The problem with the simple Fock builder is that permutational symmetries of the Fock,
    // density, and two-electron integrals are not taken into account to reduce the cost.
    // To make the simple Fock builder efficient we must rearrange our computation.
    // The most expensive step in Fock matrix construction is the evaluation of 2-e integrals;
    // hence we must minimize the number of computed integrals by taking advantage of their permutational
    // symmetry. Due to the multiplicative and Hermitian nature of the Coulomb kernel (and realness
    // of the Gaussians) the permutational symmetry of the 2-e ints is given by the following relations:
    //
    // (12|34) = (21|34) = (12|43) = (21|43) = (34|12) = (43|12) = (34|21) = (43|21)
    //
    // (here we use chemists' notation for the integrals, i.e in (ab|cd) a and b correspond to
    // electron 1, and c and d -- to electron 2).
    //
    // It is easy to verify that the following set of nested loops produces a permutationally-unique
    // set of integrals:
    // foreach a = 0 .. n-1
    //   foreach b = 0 .. a
    //     foreach c = 0 .. a
    //       foreach d = 0 .. (a == c ? b : c)
    //         compute (ab|cd)
    //
    // The only complication is that we must compute integrals over shells. But it's not that complicated ...
    //
    // The real trick is figuring out to which matrix elements of the Fock matrix each permutationally-unique
    // (ab|cd) contributes. STOP READING and try to figure it out yourself. (to check your answer see below)

    // loop over permutationally-unique set of shells
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

                    // ANSWER
                    // 1) each shell set of integrals contributes up to 6 shell sets of the Fock matrix:
                    //    F(a,b) += (ab|cd) * D(c,d)
                    //    F(c,d) += (ab|cd) * D(a,b)
                    //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
                    //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
                    //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
                    //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
                    // 2) each permutationally-unique integral (shell set) must be scaled by its degeneracy,
                    //    i.e. the number of the integrals/sets equivalent to it
                    // 3) the end result must be symmetrized
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

Matrix build_uhf_fock(const libint2::BasisSet &obs, const Matrix &D, const Matrix &Ds) {

    using libint2::Engine;
    using libint2::Operator;
    using libint2::Shell;

    const auto n = nbasis(obs.shells());
    Matrix G = Matrix::Zero(n, n);

    // construct the 2-electron repulsion integrals engine
    Engine engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 0);

    auto shell2bf = obs.shell2bf();

    const auto &buf = engine.results();


    // loop over permutationally-unique set of shells
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


real_t rhf_energy(const Matrix &D, const Matrix &H, const Matrix &F) {
    real_t energy = 0.0;
    for (auto i = 0; i < D.rows(); ++i)
        for (auto j = 0; j < D.rows(); ++j)
            energy += D(i, j) * (H(i, j) + F(i, j));
    return energy;
}

real_t uhf_energy(const Matrix &D, const Matrix &Dalpha, const Matrix &Dbeta, const Matrix &H, const Matrix &Falpha,
                  const Matrix &Fbeta) {
    real_t energy = 0.0;
    for (auto i = 0; i < D.rows(); ++i)
        for (auto j = 0; j < D.rows(); ++j)
            energy += D(i, j) * H(i, j) + Dalpha(i, j) * Falpha(i, j) + Dbeta(i, j) * Fbeta(i, j);
    return 0.5 * energy;
}

scf_results
RHF(const std::vector<libint2::Atom> &atoms, const libint2::BasisSet &obs, real_t nelectron, params config) {
    std::cout << std::endl
              << "Starting RHF calculation" << std::endl;
    scf_results results;
    auto enuc = compute_enuc(atoms);
    //    std::cout << "Nuclear repulsion energy = " << enuc << " Eh " << std::endl;
    auto ndocc = nelectron / 2;
    results.nao = nbasis(obs.shells());

    // Occupied and Virtual Orbitals
    results.no = 2 * (nelectron / 2);
    results.nv = (2 * results.nao) - results.no;
    std::cout << std::endl
              << "Number of occupied orbitals: " << results.no << std::endl
              << "Number of virtual orbitals: " << results.nv << std::endl;

    // Initializing Libint
    libint2::initialize();

    // Overlap Integrals
    auto S = compute_1body_ints(obs, libint2::Operator::overlap);
    //    std::cout << "\nOverlap Integrals:\n";
    //    std::cout << S << std::endl;

    // Kinetic Energy Integrals
    auto T = compute_1body_ints(obs, libint2::Operator::kinetic);
    //    std::cout << "\nKinetic-Energy Integrals:\n";
    //    std::cout << T << std::endl;

    // Nuclear Attraction Integrals
    Matrix V = compute_1body_ints(obs, libint2::Operator::nuclear, atoms);
    //    std::cout << "\nNuclear Attraction Integrals:\n";
    //    std::cout << V << std::endl;

    // Core Hamiltonian = T + V
    Matrix H = T + V;
    //    std::cout << "\nCore Hamiltonian:\n";
    //    std::cout << H << std::endl;

    // T and V no longer needed, free up the memory
    T.resize(0, 0);
    V.resize(0, 0);

    Matrix D;
    Matrix D_minbs = compute_soad(atoms);
    if (config.basis == "STO-3G") {
        std::cout << std::endl
                  << "Using SAD for initial guess" << std::endl;
        D = D_minbs;
    } else {
        //D = H;
        std::cout << std::endl
                  << "Building initial guess" << std::endl;
        D = density_guess(ndocc, results.nao);
    }

    //    std::cout << "\nInitial Density Matrix:\n";
    //    std::cout << D << std::endl;

    // SCF Loop
    real_t rmsd;
    real_t ediff;
    real_t ehf;

    for (auto iter = 1; iter < config.maxiter; ++iter) {
        // Save a copy of the energy and the density
        auto ehf_last = ehf;
        auto D_last = D;

        // New Fock matrix
        auto F = H;
        //F += compute_2body_fock_simple(shells, D);
        F += build_fock(obs.shells(), D);

        //        if (iter == 1) {
        //            std::cout << "\nFock Matrix:\n";
        //            std::cout << F << std::endl;
        //        }

        // solve F C = e S C
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
        auto eps = gen_eig_solver.eigenvalues();
        auto C = gen_eig_solver.eigenvectors();

        // compute density, D = C(occ) . C(occ)T
        auto C_occ = C.leftCols(ndocc);
        D = C_occ * C_occ.transpose();

        // compute HF energy
        ehf = rhf_energy(D, H, F);

        // compute difference with last iteration
        ediff = ehf - ehf_last;
        rmsd = (D - D_last).norm();


        if (iter == 1)
            std::cout << "\n\n Iter        E(elec)              E(tot)               Delta(E)             RMS(D)\n";
        printf(" %02d %20.12f %20.12f %20.12f %20.12f\n", iter, ehf, ehf + enuc, ediff, rmsd);

        if (fabs(ediff) < config.scf_conv && fabs(rmsd) < config.scf_conv) {
            results.energy = ehf + enuc;
            results.D = D;
            results.F = F;
            results.C = C;
            results.moes = eps;
            break;
        } else
            continue;
    }
    std::cout << std::endl
              << "Hartree-Fock Energy = " << results.energy << " Eh" << std::endl;

    libint2::finalize();// done with libint
    return results;
}

scf_results
UHF(const std::vector<libint2::Atom> &atoms, const libint2::BasisSet &obs, real_t nelectron, params config) {
    scf_results results;
    results.nbeta = (nelectron - config.multiplicity + 1) / 2;
    results.nalpha = results.nbeta + config.multiplicity - 1;

    std::cout << std::endl
              << "Number of alpha electrons: " << results.nalpha << std::endl
              << "Number of beta electrons: " << results.nbeta << std::endl;

    results.nao = nbasis(obs.shells());
    // Occupied and Virtual Orbitals
    results.no = 2 * (nelectron / 2);
    results.nv = (2 * results.nao) - results.no;
    std::cout << std::endl
              << "Number of occupied orbitals: " << results.no << std::endl
              << "Number of virtual orbitals: " << results.nv << std::endl;

    std::cout << std::endl
              << "Starting UHF calculation" << std::endl;
    auto enuc = compute_enuc(atoms);
    //    std::cout << "Nuclear repulsion energy = " << enuc << " Eh " << std::endl;

    // Initializing Libint
    libint2::initialize();

    // Overlap Integrals
    auto S = compute_1body_ints(obs, libint2::Operator::overlap);
    //    std::cout << "\nOverlap Integrals:\n";
    //    std::cout << S << std::endl;

    // Kinetic Energy Integrals
    auto T = compute_1body_ints(obs, libint2::Operator::kinetic);
    //    std::cout << "\nKinetic-Energy Integrals:\n";
    //    std::cout << T << std::endl;

    // Nuclear Attraction Integrals
    Matrix V = compute_1body_ints(obs, libint2::Operator::nuclear, atoms);
    //    std::cout << "\nNuclear Attraction Integrals:\n";
    //    std::cout << V << std::endl;

    // Core Hamiltonian = T + V
    Matrix H = T + V;
    //    std::cout << "\nCore Hamiltonian:\n";
    //    std::cout << H << std::endl;

    // T and V no longer needed, free up the memory
    T.resize(0, 0);
    V.resize(0, 0);

    // Building Initial Densities
    Matrix Dalpha = density_guess(results.nalpha, results.nao);
    Matrix Dbeta = density_guess(results.nbeta, results.nao);
    Matrix D = Dalpha + Dbeta;// Total Density Matrix


    // SCF Loop
    real_t rmsd;
    real_t ediff;
    real_t euhf;

    for (auto iter = 1; iter < config.maxiter; ++iter) {
        // Save copy of energy and density
        auto euhf_last = euhf;
        auto D_last = D;

        // New Fock Matrices
        auto Falpha = H;
        Falpha += build_uhf_fock(obs.shells(), D, Dalpha);
        auto Fbeta = H;
        Fbeta += build_uhf_fock(obs.shells(), D, Dbeta);

        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver1(Falpha, S);
        auto eps_alpha = solver1.eigenvalues();
        auto C_alpha = solver1.eigenvectors();

        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver2(Fbeta, S);
        auto eps_beta = solver2.eigenvalues();
        auto C_beta = solver2.eigenvectors();

        // Density Matrices
        auto Ca_occ = C_alpha.leftCols(results.nalpha);
        Dalpha = Ca_occ * Ca_occ.transpose();

        auto Cb_occ = C_beta.leftCols(results.nbeta);
        Dbeta = Cb_occ * Cb_occ.transpose();

        D = Dalpha + Dbeta;

        // UHF Energy
        euhf = uhf_energy(D, Dalpha, Dbeta, H, Falpha, Fbeta);

        // compute difference with last iteration
        ediff = euhf - euhf_last;
        rmsd = (D - D_last).norm();

        if (iter == 1)
            std::cout << "\n\n Iter        E(elec)              E(tot)               Delta(E)             RMS(D)\n";
        printf(" %02d %20.12f %20.12f %20.12f %20.12f\n", iter, euhf, euhf + enuc, ediff, rmsd);

        if (fabs(ediff) < config.scf_conv && fabs(rmsd) < config.scf_conv) {
            results.energy = euhf + enuc;
            results.D = D;
            results.Fa = Falpha;
            results.Fb = Fbeta;
            results.Ca = C_alpha;
            results.Cb = C_beta;
            results.moes_a = eps_alpha;
            results.moes_b = eps_beta;
            break;
        } else
            continue;
    }
    std::cout << std::endl
              << "Hartree-Fock Energy = " << results.energy << " Eh" << std::endl;

    return results;
}

// EOF
