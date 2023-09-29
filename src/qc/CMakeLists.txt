add_library(QC2_lib scf.cpp scf.h core.cpp core.h integrals.cpp integrals.h mppt.cpp mppt.h)
target_link_libraries(QC2_lib Libint2::cxx BTAS nlohmann_json::nlohmann_json)
