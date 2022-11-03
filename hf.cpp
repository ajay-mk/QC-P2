// Contains functions relevant to Hartree Fock Algorithm
//
// Created by Ajay Melekamburath on 10/26/22.
//

#include <iostream>
#include <string>
#include <vector>
#include <istream>

#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

// Libint Gaussian integrals library
#include <libint2.hpp>
#include <libint2/chemistry/sto3g_atomic_density.h>
#if !LIBINT2_CONSTEXPR_STATICS
#  include <libint2/statics_definition.h>
#endif

//TypeDefs
using real_t = libint2::scalar_type;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

// Structs
struct params {
    std::string type;
    std::string basis;
    int maxiter;
    real_t conv;
};

// Function Definitions

params read_config(const std::string& config_file){
    std::cout << "Reading configurations from " << config_file << std::endl;
    // For now setting some sample numbers
    params config;
    std::ifstream input (config_file);
    if (input.is_open()){
        input >> config.type;
        input >> config.basis;
        input >> config.maxiter;
        input >> config.conv;
    }

    return config;
}

// Reading Geometry from input file
std::vector<libint2::Atom> read_geometry(const std::string& filename) {
    std::cout << "Reading geometry from " << filename << std::endl;
    std::ifstream is(filename);
    assert(is.good());

    // check the extension: if .xyz, assume the standard XYZ format, otherwise throw an exception
    if ( filename.rfind(".xyz") != std::string::npos)
        return libint2::read_dotxyz(is);
    else
        throw std::invalid_argument("Only .xyz files are accepted as input");
}
// Printing Coordinates
void print_geometry(const std::vector<libint2::Atom>& atoms){
    std::cout << std::endl
              << "Molecular Geometry" << std::endl;
    for (auto i =0; i < atoms.size(); i++){
        std::cout << atoms[i].atomic_number << " " << atoms[i].x << " " << atoms[i].y << " " << atoms[i].z
                  <<std::endl;
    }
    std::cout << std::endl;
}

// Computing Nuclear Repulsion Energy
double compute_enuc(const std::vector<libint2::Atom>& atoms){
    auto num = 0.0;
    for (auto i = 0; i < atoms.size(); i++)
        for (auto j = i + 1; j < atoms.size(); j++) {
            auto xij = atoms[i].x - atoms[j].x;
            auto yij = atoms[i].y - atoms[j].y;
            auto zij = atoms[i].z - atoms[j].z;
            auto r2 = xij*xij + yij*yij + zij*zij;
            auto r = sqrt(r2);
            num += atoms[i].atomic_number * atoms[j].atomic_number / r;
        }
    return num;
}
// Counting number of basis functions
size_t nbasis(const std::vector<libint2::Shell>& shells) {
    size_t n = 0;
    for (const auto& shell: shells)
        n += shell.size();
    return n;
}

std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells) {
    std::vector<size_t> result;
    result.reserve(shells.size());

    size_t n = 0;
    for (auto shell: shells) {
        result.push_back(n);
        n += shell.size();
    }

    return result;
}

Matrix compute_1body_ints(const std::vector<libint2::Shell>& shells,
                          libint2::Operator obtype,
                          const std::vector<libint2::Atom>& atoms)
{
    using libint2::Shell;
    using libint2::Engine;
    using libint2::Operator;

    const auto n = nbasis(shells);
    Matrix result(n,n);

    // construct the overlap integrals engine
    Engine engine(obtype, max_nprim(shells), max_l(shells), 0);
    // nuclear attraction ints engine needs to know where the charges sit ...
    // the nuclei are charges in this case; in QM/MM there will also be classical charges
    if (obtype == Operator::nuclear) {
        std::vector<std::pair<real_t,std::array<real_t,3>>> q;
        for(const auto& atom : atoms) {
            q.push_back( {static_cast<real_t>(atom.atomic_number), {{atom.x, atom.y, atom.z}}} );
        }
        engine.set_params(q);
    }

    auto shell2bf = map_shell_to_basis_function(shells);

    // buf[0] points to the target shell set after every call  to engine.compute()
    const auto& buf = engine.results();

    // loop over unique shell pairs, {s1,s2} such that s1 >= s2
    // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
    for(auto s1=0; s1!=shells.size(); ++s1) {

        auto bf1 = shell2bf[s1]; // first basis function in this shell
        auto n1 = shells[s1].size();

        for(auto s2=0; s2<=s1; ++s2) {

            auto bf2 = shell2bf[s2];
            auto n2 = shells[s2].size();

            // compute shell pair
            engine.compute(shells[s1], shells[s2]);

            // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
            Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
            result.block(bf1, bf2, n1, n2) = buf_mat;
            if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
                result.block(bf2, bf1, n2, n1) = buf_mat.transpose();

        }
    }

    return result;
}
//Computes Superposition-Of-Atomic-Densities guess for the molecular density matrix
//in minimal basis; occupies subshells by smearing electrons evenly over the orbitals
Matrix compute_soad(const std::vector<libint2::Atom>& atoms) {
    // compute number of atomic orbitals
    size_t nao = 0;
    for (const auto& atom : atoms) {
        const auto Z = atom.atomic_number;
        nao += libint2::sto3g_num_ao(Z);
    }

    // compute the minimal basis density
    Matrix D = Matrix::Zero(nao, nao);
    size_t ao_offset = 0;  // first AO of this atom
    for (const auto& atom : atoms) {
        const auto Z = atom.atomic_number;
        const auto& occvec = libint2::sto3g_ao_occupation_vector(Z);
        for(const auto& occ: occvec) {
            D(ao_offset, ao_offset) = occ;
            ++ao_offset;
        }
    }

    return D * 0.5;  // we use densities normalized to # of electrons/2
}

// Fock Builder

Matrix compute_2body_fock_simple(const std::vector<libint2::Shell> &shells,
                                 const Matrix &D) {

    using libint2::Engine;
    using libint2::Operator;
    using libint2::Shell;

    const auto n = nbasis(shells);
    Matrix G = Matrix::Zero(n, n);

    // construct the electron repulsion integrals engine
    Engine engine(Operator::coulomb, max_nprim(shells), max_l(shells), 0);

    auto shell2bf = map_shell_to_basis_function(shells);

    // buf[0] points to the target shell set after every call  to engine.compute()
    const auto &buf = engine.results();

    // loop over shell pairs of the Fock matrix, {s1,s2}
    // Fock matrix is symmetric, but skipping it here for simplicity (see compute_2body_fock)
    for (auto s1 = 0; s1 != shells.size(); ++s1) {

        auto bf1_first = shell2bf[s1];// first basis function in this shell
        auto n1 = shells[s1].size();

        for (auto s2 = 0; s2 != shells.size(); ++s2) {

            auto bf2_first = shell2bf[s2];
            auto n2 = shells[s2].size();

            // loop over shell pairs of the density matrix, {s3,s4}
            // again symmetry is not used for simplicity
            for (auto s3 = 0; s3 != shells.size(); ++s3) {

                auto bf3_first = shell2bf[s3];
                auto n3 = shells[s3].size();

                for (auto s4 = 0; s4 != shells.size(); ++s4) {

                    auto bf4_first = shell2bf[s4];
                    auto n4 = shells[s4].size();

                    // Coulomb contribution to the Fock matrix is from {s1,s2,s3,s4} integrals
                    engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
                    const auto *buf_1234 = buf[0];
                    if (buf_1234 == nullptr)
                        continue;// if all integrals screened out, skip to next quartet

                    // we don't have an analog of Eigen for tensors (yet ... see github.com/BTAS/BTAS, under development)
                    // hence some manual labor here:
                    // 1) loop over every integral in the shell set (= nested loops over basis functions in each shell)
                    // and 2) add contribution from each integral
                    for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
                        const auto bf1 = f1 + bf1_first;
                        for (auto f2 = 0; f2 != n2; ++f2) {
                            const auto bf2 = f2 + bf2_first;
                            for (auto f3 = 0; f3 != n3; ++f3) {
                                const auto bf3 = f3 + bf3_first;
                                for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                                    const auto bf4 = f4 + bf4_first;
                                    G(bf1, bf2) += D(bf3, bf4) * 2.0 * buf_1234[f1234];
                                }
                            }
                        }
                    }

                    // exchange contribution to the Fock matrix is from {s1,s3,s2,s4} integrals
                    engine.compute(shells[s1], shells[s3], shells[s2], shells[s4]);
                    const auto *buf_1324 = buf[0];

                    for (auto f1 = 0, f1324 = 0; f1 != n1; ++f1) {
                        const auto bf1 = f1 + bf1_first;
                        for (auto f3 = 0; f3 != n3; ++f3) {
                            const auto bf3 = f3 + bf3_first;
                            for (auto f2 = 0; f2 != n2; ++f2) {
                                const auto bf2 = f2 + bf2_first;
                                for (auto f4 = 0; f4 != n4; ++f4, ++f1324) {
                                    const auto bf4 = f4 + bf4_first;
                                    G(bf1, bf2) -= D(bf3, bf4) * buf_1324[f1324];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return G;
}

Matrix compute_2body_fock(const std::vector<libint2::Shell> &shells,
                          const Matrix &D) {

    using libint2::Engine;
    using libint2::Operator;
    using libint2::Shell;

    auto time_elapsed = std::chrono::duration<double>::zero();

    const auto n = nbasis(shells);
    Matrix G = Matrix::Zero(n, n);

    // construct the 2-electron repulsion integrals engine
    Engine engine(Operator::coulomb, max_nprim(shells), max_l(shells), 0);

    auto shell2bf = map_shell_to_basis_function(shells);

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
    for (auto s1 = 0; s1 != shells.size(); ++s1) {

        auto bf1_first = shell2bf[s1];// first basis function in this shell
        auto n1 = shells[s1].size();  // number of basis functions in this shell

        for (auto s2 = 0; s2 <= s1; ++s2) {

            auto bf2_first = shell2bf[s2];
            auto n2 = shells[s2].size();

            for (auto s3 = 0; s3 <= s1; ++s3) {

                auto bf3_first = shell2bf[s3];
                auto n3 = shells[s3].size();

                const auto s4_max = (s1 == s3) ? s2 : s3;
                for (auto s4 = 0; s4 <= s4_max; ++s4) {

                    auto bf4_first = shell2bf[s4];
                    auto n4 = shells[s4].size();

                    // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
                    auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
                    auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
                    auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
                    auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

                    const auto tstart = std::chrono::high_resolution_clock::now();

                    engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
                    const auto *buf_1234 = buf[0];
                    if (buf_1234 == nullptr)
                        continue;// if all integrals screened out, skip to next quartet

                    const auto tstop = std::chrono::high_resolution_clock::now();
                    time_elapsed += tstop - tstart;

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
