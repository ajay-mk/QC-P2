//
// Created by Ajay Melekamburath on 11/24/22.
//

#ifndef P2_GENERAL_H
#define P2_GENERAL_H

// Some general functions

#include <iostream>
#include <vector>

// Libint Gaussian integrals library
#include <libint2.hpp>
#if !LIBINT2_CONSTEXPR_STATICS
#include <libint2/statics_definition.h>
#endif

// Typedefs
using real_t = libint2::scalar_type;


// Structs
struct params {
    std::string inputfile;
    std::string type;
    std::string scf;
    std::string basis;
    double multiplicity;
    int maxiter;
    real_t conv;
};

std::vector<libint2::Atom> read_geometry(const std::string &filename);
void print_geometry(const std::vector<libint2::Atom> &atoms);
params read_config(const std::string& config_file);
size_t nbasis(const std::vector<libint2::Shell>& shells);

#endif//P2_GENERAL_H

// EOF