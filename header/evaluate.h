#pragma once
#include "seal/seal.h"
#include "CKKS_params.h"
#include "parse.h"
#include <iomanip>
#include <cmath>
using namespace std;
using namespace seal;

void evaluate_helut(CKKS_params& pms, HELUT_Info& hi, double e, double p, bool mid_print, bool res_print);
void evaluate_multi_remez(CKKS_params& pms, MR_Info& entries, double e, double p, bool mid_print, bool res_print);

Ciphertext evaluate_func(CKKS_params& pms, vector<double> coeff, Ciphertext& x);
// Ciphertext evaluate_func_even(CKKS_params& pms, vector<double> coeff, double scale, Ciphertext& x, bool correct_factor);
Ciphertext eval_even(CKKS_params& pms, vector<double> coeff, vector<double> correction_factors, Ciphertext& x);
Ciphertext eval_odd(CKKS_params& pms, vector<double> coeff, vector<double> correction_factors, Ciphertext& x);
vector<double> calculate_correction_factor(CKKS_params& pms, vector<double> coeff, Ciphertext& x);
vector<double> calculate_correction_factor_odd(CKKS_params& pms, vector<double> coeff, Ciphertext& x);
double calculate_cleanse_correction_factor(CKKS_params& pms, Ciphertext& x);
Ciphertext cleanse(CKKS_params& pms, Ciphertext& x);
int check_coeff_type(vector<double> coeff);
int evaluate_depth(vector<double> coeff, bool include_plain);

std::vector<double> linspace(double a, double b, std::size_t n);
void compare_result(vector<double> res1, vector<double> res2, string title1, string title2, int size);
double compare_result_log(const std::vector<double>& res, int size, bool print);


// Plain data 전용
std::vector<double> exp_vector(const std::vector<double>& x, int exp);
std::vector<double> mult_vector(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> add_vector(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> mult_vector_plain(const std::vector<double>& a, double scalar);
std::vector<double> add_vector_plain(const std::vector<double>& a, double scalar);
std::vector<double> evaluate_function(const std::vector<double>& coeff, const std::vector<double>& x);
