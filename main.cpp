#include <iostream>
#include <vector>
#include <fstream>


// Libint Gaussian integrals library
#include <libint2.hpp>
#if !LIBINT2_CONSTEXPR_STATICS
#  include <libint2/statics_definition.h>
#endif

using std::cout, std::endl;

// Functions
std::vector<libint2::Atom> read_geometry(const std::string& filename);


int main(int argc, char* argv[]) {
    using std::cout, std::endl, std::cerr;

    // Input Geometry
    const auto filename = argv[1];
    // Calculation Configurations
    //const auto config = argv[2];

    std::vector<libint2::Atom> atoms = read_geometry(filename);

    return 0;
}

// Function Definitions
std::vector<libint2::Atom> read_geometry(const std::string & filename){
    cout << "Reading geometry from " << filename << endl;
    std::ifstream input(filename);
    if(filename.rfind(".xyz"))
        return libint2::read_dotxyz(input);
    else
        throw std::invalid_argument("Only .xyz files accepted as input");
}


