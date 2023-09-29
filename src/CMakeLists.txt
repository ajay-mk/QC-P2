add_library(qc scf.cpp scf.h core.cpp core.h integrals.cpp integrals.h mppt.cpp mppt.h)
target_link_libraries(qc Libint2::cxx BTAS nlohmann_json::nlohmann_json)
#target_compile_features(qc PUBLIC cxx_std_11)
