#include <iomanip>
#include <iostream>
#include <vector>
#include "../src/utils.h"

int main(int argc, char *argv[]) {
  using namespace qc;

  std::cout << std::setprecision(12);  // prints 12 decimals

  const auto config = utils::read_config(argv[1]);
  const auto geom = utils::read_geometry(config.geom_file);
}
