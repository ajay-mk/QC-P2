//
// Created by Ajay Melekamburath on 11/14/22.
//
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
    std::string type;
    std::string basis;
    double multiplicity;
    int maxiter;
    real_t conv;
};

// Function Definitions

params read_config(const std::string& config_file){
    std::cout << "Reading configurations from " << config_file << std::endl;
    params config;
    // Expected Format of Config File
    // Method
    // Basis Set
    // Multiplicity
    // SCF Max. Iter.
    // SCF Conv
    std::ifstream input (config_file);
    if (input.is_open()){
        input >> config.type;
        input >> config.basis;
        input >> config.multiplicity;
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


// EOF