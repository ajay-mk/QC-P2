//
// Created by Ajay Melekamburath on 5/13/23.
//

#include "core.h"

namespace qc::core {
std::vector<libint2::Atom> read_geometry(const std::string &filename) {
  std::cout << "Reading geometry from " << filename << std::endl;
  std::ifstream is(filename);
  assert(is.good() && "Bad input file.");

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

template <typename T>
T read_json_item(const nlohmann::json &json, const std::string &key,
                 const T &def_value) {
  if (json.contains(key)) return json[key];
  return def_value;
}

input read_config(const std::string &config_file) {
  using nlohmann::json;
  // TODO: Switch to Boost::json
  std::cout << "\nReading configuration from " << config_file << std::endl;
  input config;  // make the struct
  std::ifstream f(config_file);
  json input = json::parse(f);

  // new implementation
  config.geom_file = read_json_item<std::string>(input, "filename", "h2o.xyz");
  config.type = read_json_item<std::string>(input, "type", "RHF");
  config.basis = read_json_item<std::string>(input, "basis", "STO-3G");
  config.maxiter = read_json_item<int>(input, "maxiter", 50);
  config.scf_conv = read_json_item<real_t>(input, "scf_conv", 1e-8);
  config.cc_conv = read_json_item<real_t>(input, "cc_conv", config.scf_conv);
  config.charge = read_json_item<int>(input, "charge", 0);
  config.multiplicity = read_json_item<int>(input, "multiplicity", 1);

  // SCF
  config.ref = read_json_item<std::string>(input, "ref", "RHF");
  assert(config.ref == "RHF" || config.ref == "rhf" || config.ref == "UHF" ||
         config.ref == "uhf" && "Unsupported reference in input");

  // try and catch block for cc_conv <= scf_conv
  try {
    if (config.cc_conv <= config.scf_conv) {
      throw std::invalid_argument(
          "Warning: requested cc_conv is higher than scf_conv. "
          "Setting cc_conv = scf_conv.");
    }
  } catch (const std::invalid_argument &e) {
    std::cerr << e.what() << std::endl;
    config.cc_conv = config.scf_conv;
  }
  std::cout << "Input JSON:\n" << input.dump(5) << std::endl;
  return config;
}

}  // namespace qc::core

namespace qc::utils {

size_t nbasis(const std::vector<libint2::Shell> &shells) {
  size_t n = 0;
  for (const auto &shell : shells) n += shell.size();
  return n;
}
}  // namespace qc::utils
