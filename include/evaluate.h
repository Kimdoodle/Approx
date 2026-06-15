#pragma once
#include "seal/seal.h"
#include "CKKS_params.h"
#include "polyEval_class.h"
#include "parse.h"
#include "error_bound.h"
#include "fakeciphertext.h"
#include "eval_plain.h"
#include <math.h>
#include <iomanip>

// using CT = std::tuple<std::vector<double>, seal::Ciphertext, fct::Ciphertext>;


// Core - HELUT/REMEZ evaluation
void evaluate_helut(CKKS_params& pms, vector<double>& x, bool test_time);
void evaluate_multi_remez(CKKS_params& pms, vector<double>& x, std::vector<std::shared_ptr<Decomp>> dcmps, bool test_time);

// Submodules
double_t correction_factor_HELUT(CKKS_params& pms, seal::Ciphertext& x);
double_t correction_factor_REMEZ(CKKS_params& pms, std::vector<double>& coeff, seal::Ciphertext& ct);
seal::Ciphertext evaluate_xi(CKKS_params& pms, seal::Ciphertext& ct, XI xi, std::map<int, seal::Ciphertext> powers, double_t correction_factor);
seal::Ciphertext evaluate_polynomial(CKKS_params& pms, seal::Ciphertext& ct, std::shared_ptr<Decomp> dcmp, std::map<int, seal::Ciphertext> powers, double_t correction_factor=1.0);

// Evaluation modules - PT, CT, FCT
// CT eval_poly_ct(CKKS_params& pms, shared_ptr<Decomp> dcmp, CT x, ErrBound& eb, int copy_count, bool step_debug=false, bool res_debug=true);
// seal::Ciphertext eval_poly_ct(CKKS_params& pms, shared_ptr<Decomp> dcmp, seal::Ciphertext& x);
// CT eval_poly_cl(CKKS_params& pms, CT x, ErrBound& eb, int copy_count, bool step_debug=false, bool res_debug=true);

namespace ct
{
    seal::Ciphertext evaluate(seal::Ciphertext& x, EvalStep es, CKKS_params& pms);
}
// Comparison modules
void compare_boundary(std::vector<double> ptx, seal::Ciphertext ctx, std::pair<double, double> fct_interval, CKKS_params& pms, int copy_count);
double compare_precision(std::vector<double> ptx, seal::Ciphertext ctx, std::pair<double, double> fct_interval, CKKS_params& pms, int copy_count);

// void compare_result_pt(const std::vector<double> &ct_res, const std::vector<double> &pt_res, int size);
// void compare_result_IF(const std::vector<double>& res, int size);


// std::vector<double> evaluate_polynomial_dcmp(const std::vector<double>& coeff, const std::vector<double>& x, EvalStep es);
