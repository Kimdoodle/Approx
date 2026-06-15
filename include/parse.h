#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <stdexcept>
#include <map>
#include <limits>

#include "json.hpp"
#include "polyEval_class.h"
#include "error_bound.h"

std::vector<std::vector<double>> parse_remez_coeff(int p_num, int e_num, std::string critieria, bool print = false);
int calcuate_remez_depth(std::vector<std::vector<double>> coeffs);

nlohmann::json load_decomp_cache(const std::string& filename);

// Return coeff key.
std::string get_poly_type_key(const std::vector<double>& coeff);

// restore decomp class using json file.
std::shared_ptr<Decomp> reconstruct_decomp_from_cache(const json& plan, Poly poly, int depth=0);

// std::variant<seal::Plaintext, seal::Ciphertext> get_value2(std::string key, std::map<int, seal::Ciphertext>& powers, std::map<int, seal::Plaintext>& coeffs, std::map<int, seal::Ciphertext>& terms);
// std::variant<seal::Plaintext, seal::Ciphertext> operate(std::string key1, std::string key2, char op, std::map<int, seal::Ciphertext>& powers, std::map<int, seal::Plaintext>& coeffs, std::map<int, seal::Ciphertext>& terms);