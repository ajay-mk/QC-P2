#include <iomanip>
#include <iostream>

#include <qc/core.h>"
#include <qc/mppt.h>
#include <qc/scf.h>

int main(int argc, char *argv[]) {
  using namespace qc;

  std::cout << std::setprecision(12);  // prints 12 decimals

  const auto config = core::read_config(argv[1]);
  const auto atoms = core::read_geometry(config.geom_file);

  //  auto scf = SCF(atoms, config);
  auto mp2 = MP2(atoms, config);
}
