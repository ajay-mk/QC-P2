//
// Created by Ajay Melekamburath on 5/13/23.
//

#include "utils.h"

namespace qc::utils {
std::vector<libint2::Atom> read_geometry(const std::string &filename) {
  std::cout << "\nReading geometry from " << filename << std::endl;
  std::ifstream is(filename);
  assert(is.good() && "Bad input file.");

  // check the extension: if .xyz, assume the standard XYZ format, otherwise
  // throw an exception
  if (filename.rfind(".xyz") != std::string::npos)
    return libint2::read_dotxyz(is);
  else
    throw std::invalid_argument("Only .xyz files are accepted as input");
}

void print_geometry(const std::vector<libint2::Atom> &atoms) {
  std::cout << "\nGeometry: " << std::endl;
  for (const auto &atom : atoms)
    std::cout << atom.atomic_number << " " << atom.x << " " << atom.y << " "
              << atom.z << std::endl;
}

size_t nbasis(const std::vector<libint2::Shell> &shells) {
  size_t n = 0;
  for (const auto &shell : shells) n += shell.size();
  return n;
}

template <typename T>
T read_json_item(const nlohmann::json &json, const std::string &key,
                 const T &def_value) {
  if (json.contains(key)) return json[key];
  return def_value;
}

input read_config(const std::string &config_file) {
  using nlohmann::json;
  std::cout << "\nReading configuration from " << config_file << std::endl;
  input config;  // make the struct
  std::ifstream f(config_file);
  json input = json::parse(f);

  // new implementation
  config.geom_file = read_json_item<std::string>(input, "filename", "h2o.xyz");
  config.type = read_json_item<std::string>(input, "type", "RHF");
  config.basis = read_json_item<std::string>(input, "basis", "3-21G");
  config.maxiter = read_json_item<int>(input, "maxiter", 50);
  config.scf_conv = read_json_item<double>(input, "scf_conv", 1e-8);
  config.cc_conv = read_json_item<double>(input, "cc_conv", config.scf_conv);
  config.charge = read_json_item<int>(input, "charge", 0);
  config.multiplicity = read_json_item<int>(input, "multiplicity", 1);

  // SCF
  config.ref = read_json_item<std::string>(input, "ref", "RHF");
  assert(config.ref == "RHF" || config.ref == "rhf" || config.ref == "UHF" ||
         config.ref == "uhf" && "Unsupported reference in input");

  assert(config.cc_conv >= config.scf_conv &&
         "Requested CC convergence is higher than SCF convergence");

  std::cout << input.dump(4) << std::endl;
  return config;
}

}  // namespace qc::utils
