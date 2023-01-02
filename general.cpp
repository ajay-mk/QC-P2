//
// Created by Ajay Melekamburath on 11/14/22.
//
// Some general functions

#include <iostream>

#include "general.h"

// Function Definitions
///TODO: Dictionary type construct for config
params read_config(const std::string& config_file){
    std::cout << std::endl
              << "Reading configuration from " << config_file << std::endl;
    params config;
    // Expected Format of Config File
    // Input Geometry
    // Method
    // Basis Set
    // Multiplicity
    // SCF Max. Iter.
    // SCF Conv
    std::ifstream input (config_file);
    if (input.is_open()){
        input >> config.inputfile;
        input >> config.type;
        input >> config.ref;
        input >> config.basis;
        input >> config.multiplicity;
        input >> config.maxiter;
        input >> config.conv;
    }

    return config;
}

// Reading Geometry from input file
std::vector<libint2::Atom> read_geometry(const std::string& filename) {
    std::cout <<
            std::endl << "Reading geometry from " << filename << std::endl;
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
              << "Geometry: " << std::endl;
    for (auto i =0; i < atoms.size(); i++){
        std::cout << atoms[i].atomic_number << " " << atoms[i].x << " " << atoms[i].y << " " << atoms[i].z
                  <<std::endl;
    }
    std::cout << std::endl;
}

// Counting number of basis functions
size_t nbasis(const std::vector<libint2::Shell>& shells) {
    size_t n = 0;
    for (const auto& shell: shells)
        n += shell.size();
    return n;
}

// EOF