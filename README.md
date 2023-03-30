### Quantum Chemistry Codes using C++

#### Prerequisites

- [Libint](https://github.com/evaleev/libint)
- [BTAS](https://github.com/ValeevGroup/BTAS)
- [nlohmann/json](https://github.com/nlohmann/json)

All dependencies except Libint will be fetched and compiled by CMake.

Install Libint and add to `-DCMAKE_PREFIX_PATH` while building.

#### Features

- Restricted Hartree-Fock
- Unrestricted Hartree-Fock
- Second Order Møller–Plesset Perturbation Theory (MP2)
- Coupled Cluster Singles and Doubles, Perturbative Triples (CCSD & CCSD(T))
