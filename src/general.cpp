//
// Created by Ajay Melekamburath on 11/14/22.
//
// Some general functions

#include <iostream>

#include "general.h"

// Function Definitions

params read_config_json(const std::string &config_file) {
    using nlohmann::json;

    std::cout << std::endl
              << "Reading configuration from " << config_file << std::endl;
    params config;

    json input;
    std::ifstream stream(config_file);
    stream >> input;

    // Required Parameters
    //TODO: Throw error if the below parameters are empty
    config.inputfile = input["filename"];
    config.type = input["type"];
    config.basis = input["basis"];
    config.maxiter = input["maxiter"];
    config.scf_conv = input["scf_conv"];

    // Optional Parameters
    // If some parameters are empty
    if (input["ref"].empty())
        config.ref = "RHF";
    else
        config.ref = input["ref"];

    if (input["charge"].empty())
        config.charge = 0;
    else
        config.charge = input["charge"];

    if (input["multiplicity"].empty())
        config.multiplicity = 1;
    else
        config.multiplicity = input["multiplicity"];

    if (input["cc_conv"].empty())
        config.cc_conv = input["scf_conv"];
    else
        config.cc_conv = input["cc_conv"];

    return config;
}

// Reading Geometry from input file
std::vector<libint2::Atom> read_geometry(const std::string &filename) {
    std::cout << std::endl
              << "Reading geometry from " << filename << std::endl;
    std::ifstream is(filename);
    assert(is.good());

    // check the extension: if .xyz, assume the standard XYZ format, otherwise throw an exception
    if (filename.rfind(".xyz") != std::string::npos)
        return libint2::read_dotxyz(is);
    else
        throw std::invalid_argument("Only .xyz files are accepted as input");
}

// Printing Coordinates
void print_geometry(const std::vector<libint2::Atom> &atoms) {
    std::cout << std::endl
              << "Geometry: " << std::endl;
    for (auto i = 0; i < atoms.size(); ++i) {
        std::cout << atoms[i].atomic_number << " " << atoms[i].x << " " << atoms[i].y << " " << atoms[i].z
                  << std::endl;
    }
    std::cout << std::endl;
}

// Counting number of basis functions
size_t nbasis(const std::vector<libint2::Shell> &shells) {
    size_t n = 0;
    for (const auto &shell: shells)
        n += shell.size();
    return n;
}

// EOF
