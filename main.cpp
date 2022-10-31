#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigenvalues>

// Libint Gaussian integrals library
#include <libint2.hpp>
#if !LIBINT2_CONSTEXPR_STATICS
#include <libint2/statics_definition.h>
#endif

// Calculation Parameters (Will add function to read from config file)
std::string basis = "STO-3G";

// Typedefs
using real_t = libint2::scalar_type;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

// Functions
std::vector<libint2::Atom> read_geometry(const std::string &filename);
void print_geometry(const std::vector<libint2::Atom> &atoms);
size_t nbasis(const std::vector<libint2::Shell> &shells);
std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell> &shells);

Matrix compute_soad(const std::vector<libint2::Atom> &atoms);
double compute_enuc(const std::vector<libint2::Atom> &atoms);
Matrix compute_1body_ints(const std::vector<libint2::Shell> &shells,
                          libint2::Operator t,
                          const std::vector<libint2::Atom> &atoms = std::vector<libint2::Atom>());

// simple-to-read, but inefficient Fock builder; computes ~16 times as many ints as possible
Matrix compute_2body_fock_simple(const std::vector<libint2::Shell> &shells,
                                 const Matrix &D);
// an efficient Fock builder; *integral-driven* hence computes permutationally-unique ints once
Matrix compute_2body_fock(const std::vector<libint2::Shell> &shells,
                          const Matrix &D);

int main(int argc, char *argv[]) {

    using std::cerr;
    using std::cout;
    using std::endl;

    using libint2::BasisSet;
    using libint2::Engine;
    using libint2::Operator;
    using libint2::Shell;

    cout << std::setprecision(12);

    // Reading geometry from input file
    const auto filename = argv[1];
    std::vector<libint2::Atom> atoms = read_geometry(filename);

    // Counting the number of electrons
    auto nelectron = 0;
    for (auto i = 0; i < atoms.size(); ++i)
        nelectron += atoms[i].atomic_number;
    const auto ndocc = nelectron / 2;

    print_geometry(atoms);
    cout << "Number of electrons = " << nelectron << endl;


    auto enuc = compute_enuc(atoms);
    cout << "Nuclear repulsion energy = " << enuc << " Eh " << endl;

    BasisSet obs(basis, atoms);
    auto nao = nbasis(obs.shells());
    cout << endl
         << "Basis Set: " << basis << endl
         << "Number of basis functions = " << nao << endl;

    // Initializing Libint
    libint2::initialize();

    // Overlap Integrals
    auto S = compute_1body_ints(obs.shells(), Operator::overlap);
    cout << "\n\tOverlap Integrals:\n";
    cout << S << endl;

    // Kinetic Energy Integrals
    auto T = compute_1body_ints(obs.shells(), Operator::kinetic);
    cout << "\n\tKinetic-Energy Integrals:\n";
    cout << T << endl;

    // Nuclear Attraction Integrals
    Matrix V = compute_1body_ints(obs.shells(), Operator::nuclear, atoms);
    cout << "\n\tNuclear Attraction Integrals:\n";
    cout << V << endl;

    // Core Hamiltonian = T + V
    Matrix H = T + V;
    cout << "\n\tCore Hamiltonian:\n";
    cout << H << endl;

    // T and V no longer needed, free up the memory
    T.resize(0, 0);
    V.resize(0, 0);

    Matrix D;
    D = compute_soad(atoms);

    cout << "\n\tInitial Density Matrix:\n";
    cout << D << endl;

    // SCF Loop
    const auto maxiter = 100;
    const real_t conv = 1e-12;
    auto iter = 0;
    real_t rmsd = 0.0;
    real_t ediff = 0.0;
    real_t ehf = 0.0;
    do {
        ++iter;

        // Save a copy of the energy and the density
        auto ehf_last = ehf;
        auto D_last = D;

        // New Fock matrix
        auto F = H;
        //F += compute_2body_fock_simple(shells, D);
        F += compute_2body_fock(obs.shells(), D);

        if (iter == 1) {
            cout << "\n\tFock Matrix:\n";
            cout << F << endl;
        }

        // solve F C = e S C
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
        auto eps = gen_eig_solver.eigenvalues();
        auto C = gen_eig_solver.eigenvectors();

        // compute density, D = C(occ) . C(occ)T
        auto C_occ = C.leftCols(ndocc);
        D = C_occ * C_occ.transpose();

        // compute HF energy
        ehf = 0.0;
        for (auto i = 0; i < nao; i++)
            for (auto j = 0; j < nao; j++)
                ehf += D(i, j) * (H(i, j) + F(i, j));

        // compute difference with last iteration
        ediff = ehf - ehf_last;
        rmsd = (D - D_last).norm();


        if (iter == 1)
            std::cout << "\n\n Iter        E(elec)              E(tot)               Delta(E)             RMS(D)\n";
        printf(" %02d %20.12f %20.12f %20.12f %20.12f\n", iter, ehf, ehf + enuc,
               ediff, rmsd);

    } while (((fabs(ediff) > conv) || (fabs(rmsd) > conv)) && (iter < maxiter));

    cout << endl
         << "Hartree-Fock Energy = " << ehf + enuc << endl;

    libint2::finalize();// done with libint

}

