#include "header/evaluate.h"
#include <math.h>


void evaluate_helut(CKKS_params& pms, HELUT_Info& hi, double e, double p, bool mid_print, bool res_print)
{
    // cout << "r, s : " << hi.r << ", " << hi.s << endl;
    vector<double> x = {0.0, e};
    auto lin = linspace(1 - e, p-1, pms.parms->poly_modulus_degree()/2-2);
    x.insert(x.end(), lin.begin(), lin.end());
    vector<double> real_result = x;
    vector<double> ckks_result;
    Ciphertext ct = pms.encrypt(x);
    Plaintext pt;

    // Calculate correction factor
    // const auto& coeff_modulus = pms.context->get_context_data(ct.parms_id())->parms().coeff_modulus();
    // int coeff_modulus_index = coeff_modulus.size() - 1;
    // double_t correction_factor = 1.0;
    // double_t scale = ct.scale();

    // for(int i=0; i<(2 + hi.r) + (2 * hi.s); i++)
    // {
    //     correction_factor *= static_cast<double_t>(coeff_modulus[coeff_modulus_index].value()) / scale;
    //     --coeff_modulus_index;
    // }
    
    // pms.encoder->encode(1.0, ct.parms_id(), correction_factor, pt);
    // pms.eva->multiply_plain_inplace(ct, pt);

    // 1. SqMethod - (1-2x^2/p^2)^(2^r)
    {
        // Base 연산 - Level 2 
        pms.eva->square_inplace(ct);
        pms.eva->relinearize_inplace(ct, pms.rlk);
        pms.eva->rescale_to_next_inplace(ct);

        pms.encoder->encode(-2.0/(p*p), ct.parms_id(), ct.scale(), pt);
        pms.eva->multiply_plain_inplace(ct, pt);
        pms.eva->rescale_to_next_inplace(ct);

        pms.encoder->encode(1.0, ct.parms_id(), ct.scale(), pt);
        pms.eva->add_plain_inplace(ct, pt);

        //debug
        if(res_print)
        {
            real_result = exp_vector(real_result, 2);
            real_result = mult_vector_plain(real_result, (-2.0/(p*p)));
            real_result = add_vector_plain(real_result, 1.0);
            if(mid_print)
            {
                ckks_result = pms.decode_ctxt(ct);
                cout << "SqMethod base result:"<<endl;
                // compare_result(real_result, ckks_result, "REAL", "CKKS", real_result.size());
                double max_err = compare_result_log(ckks_result, x.size(), true);
                cout << "---------------------------------" << endl;
            }
        }

        // r번 제곱
        for(int i=0; i<hi.r; i++)
        {
            pms.eva->square_inplace(ct);
            pms.eva->relinearize_inplace(ct, pms.rlk);
            pms.eva->rescale_to_next_inplace(ct);

            //debug
            if(res_print)
            {
                real_result = exp_vector(real_result, 2);
                if(mid_print)
                {
                    ckks_result = pms.decode_ctxt(ct);
                    cout << "SqMethod " << i+1 <<endl;
                    // compare_result(real_result, ckks_result, "REAL", "CKKS", real_result.size());
                    double max_err = compare_result_log(ckks_result, x.size(), true);
                    cout << "---------------------------------" << endl;
                }
            }
        }
    }


    // 2. Cleanse - -2x^3 + 3x^2
    vector<double> coeff_cleanse = {0.0, 0.0, 3.0, -2.0};
    for(int i=0; i<hi.s; i++)
    {
        ct = cleanse(pms, ct);
        if(res_print)
        {
            real_result = evaluate_function(coeff_cleanse, real_result);
            if(mid_print)
            {
                ckks_result = pms.decode_ctxt(ct);
                cout << "Cleanse " << i+1 << endl;
                // compare_result(real_result, ckks_result, "REAL", "CKKS", real_result.size());
                double max_err = compare_result_log(ckks_result, x.size(), true);
                cout << "---------------------------------" << endl;
            }
        }
    }

    // pms.encoder->encode(1.0, ct.parms_id(), correction_factor, pt);
    // pms.eva->multiply_plain_inplace(ct, pt);
    // ct.scale() = scale;

    if(res_print)
    {
        cout << "Final Result: " << endl;
        ckks_result = pms.decode_ctxt(ct);
        // compare_result(real_result, ckks_result, "REAL", "CKKS", real_result.size());
        compare_result_log(ckks_result, x.size(), true);
    }
}

void evaluate_multi_remez(CKKS_params& pms, MR_Info& mi, double e, double p, bool mid_print, bool res_print)
{   
    vector<double> x;
    for(int i=0; i<p; i++)
        x.push_back(i);
    // auto lin = linspace(1 - e, p-1, pms.parms->poly_modulus_degree()/2-2);
    // x.insert(x.end(), lin.begin(), lin.end());

    Ciphertext ct = pms.encrypt(x);
    vector<vector<double>> coeffs = mi.coeffs;
    vector<double> real_result, ckks_result, coeff;

    // 1. coeff
    for(int en=0; en<mi.n - mi.s; en++)
    {
        coeff = coeffs[en];
        vector<double> corr_factor = calculate_correction_factor(pms, coeff, ct);
        ct = eval(pms, coeff, corr_factor, ct);
        if(res_print)
        {
            real_result = evaluate_function(coeff, x);
            if(mid_print)
            {
                ckks_result = pms.decode_ctxt(ct);
                cout << "Coeff Evaluation " << en+1 << ", degree " << coeff.size()-1 << endl;
                compare_result(real_result, ckks_result, "REAL", "CKKS", x.size());
                // double max_err = compare_result_log(ckks_result, x.size(), true);
                cout << "------------\n";
            }
        }
        x = real_result;
    }

    // 2. Cleanse
    vector<double> coeff_cleanse = {0.0, 0.0, 3.0, -2.0};
    for(int i=0; i<mi.s; i++)
    {
        ct = cleanse(pms, ct);
        if(res_print)
        {
            real_result = evaluate_function(coeff_cleanse, real_result);
            if(mid_print)
            {
                ckks_result = pms.decode_ctxt(ct);
                cout << "Cleanse result" << endl;
                // compare_result(real_result, ckks_result, "REAL", "CKKS", x.size());
                double max_err = compare_result_log(ckks_result, x.size(), true);
                cout << "------------\n";
            }
        }
    }

    if(res_print)
    {
        cout << "Final Result: " << endl;
        ckks_result = pms.decode_ctxt(ct);
        compare_result(real_result, ckks_result, "REAL", "CKKS", real_result.size());
        // compare_result_log(ckks_result, x.size(), true);
    }
}

// 다항식 평가
Ciphertext evaluate_func(CKKS_params& pms, vector<double> coeff, Ciphertext& x)
{
    const auto& coeff_modulus = pms.context->get_context_data(x.parms_id())->parms().coeff_modulus();
    int coeff_modulus_index = coeff_modulus.size() - 1;
    Ciphertext result;
    Plaintext pt;
    int index = coeff.size() - 1;

    // 첫 단계
    double_t correction_factor = 1.0;
    double_t scale = x.scale();

    for(size_t i=0; i<index; i++)
    {
        correction_factor *= static_cast<double_t>(coeff_modulus[coeff_modulus_index].value()) / scale;
        --coeff_modulus_index;
    }

    pms.encoder->encode(coeff[index], x.parms_id(), scale * correction_factor, pt);
    pms.eva->multiply_plain(x, pt, result);
    pms.eva->rescale_to_next_inplace(result);
    // coeff_modulus_index--;

    if(coeff[index-1] != 0.0)
    {
        pms.encoder->encode(coeff[index-1], result.parms_id(), result.scale(), pt);
        pms.eva->add_plain_inplace(result, pt);
    }

    for(int i=index-2; i>=0; --i)
    {
        pms.eva->mod_reduce_to_next_inplace(x);

        pms.eva->multiply(x, result, result);
        pms.eva->relinearize_inplace(result, pms.rlk);
        pms.eva->rescale_to_next_inplace(result);
        if(coeff[i] != 0.0)
        {
            pms.encoder->encode(coeff[i], result.parms_id(), result.scale(), pt);
            pms.eva->add_plain_inplace(result, pt);       
        }
    }
    cout<< result.scale() <<'\n';
    cout<< scale <<'\n';
    result.scale() = scale;

    return result;
}

// 다항식 평가(짝수차 함수)
// ex) ax^4 + bx^2 + c = ((ax^2 + b)x^2 + c)
Ciphertext evaluate_func_even(CKKS_params& pms, vector<double> coeff, double scale, Ciphertext& x, bool correct_factor)
{
    // for(int i=0; i<coeff.size(); i++)
    //     coeff[i] *= scale;
    const auto& coeff_modulus = pms.context->get_context_data(x.parms_id())->parms().coeff_modulus();
    int coeff_modulus_index = coeff_modulus.size() - 1;
    Ciphertext x_square, result;
    Plaintext pt;
    int index = coeff.size() - 1;

    double_t correction_factor = 1.0;
    double_t init_scale = x.scale();

    for(size_t i=0; i<index/2+1; i++)
    {
        correction_factor *= static_cast<double_t>(coeff_modulus[coeff_modulus_index].value()) / init_scale;
        coeff_modulus_index--;
    }

    if(correct_factor)
    {
        pms.encoder->encode(1.0, x.parms_id(), correct_factor, pt);
        pms.eva->multiply_plain_inplace(x, pt);
    }

    // 1. x^2 계산
    pms.eva->square(x, x_square);
    pms.eva->relinearize_inplace(x_square, pms.rlk);
    pms.eva->rescale_to_next_inplace(x_square);

    // 2. 반복
    pms.encoder->encode(coeff[index], x_square.parms_id(), 1.0, pt);
    pms.eva->multiply_plain(x_square, pt, result);
    
    if(coeff[index-2] != 0.0)
    {
        pms.encoder->encode(coeff[index-2], result.parms_id(), result.scale(), pt);
        pms.eva->add_plain_inplace(result, pt);
    }

    for(int i=index-4; i>=0; i-=2)
    {
        pms.eva->multiply(x_square, result, result);
        pms.eva->relinearize_inplace(result, pms.rlk);
        pms.eva->rescale_to_next_inplace(result);
        if(coeff[i] != 0.0)
        {
            pms.encoder->encode(coeff[i], result.parms_id(), result.scale(), pt);
            pms.eva->add_plain_inplace(result, pt);
        }
        pms.eva->mod_reduce_to_next_inplace(x_square);
    }

    // pms.encoder->encode(1.0 / scale, result.parms_id(), result.scale() * correction_factor, pt);
    // pms.eva->multiply_plain_inplace(result, pt);
    // pms.eva->rescale_to_next_inplace(result);

    if(correct_factor)
    {
        cout<< result.scale() <<'\n';
        cout<< init_scale <<'\n';
        result.scale() = init_scale;
    }

    return result;
}

Ciphertext eval(CKKS_params& pms, vector<double> coeff, vector<double> correction_factors, Ciphertext& x)
{
    int index = coeff.size() - 1;
    int cf_index = 0;
    double scale = x.scale();
    Ciphertext x_square, temp, result;
    Plaintext pt;

    //x^2
    pms.eva->square(x, x_square);
    pms.eva->relinearize_inplace(x_square, pms.rlk);
    pms.eva->rescale_to_next_inplace(x_square);

    // 첫 번째 항 = (a*f1) * x^2 + b
    pms.encoder->encode(coeff[index], x_square.parms_id(), scale * correction_factors[cf_index++], pt);
    pms.eva->multiply_plain(x_square, pt, result);
    pms.eva->rescale_to_next_inplace(result);

    pms.encoder->encode(coeff[index-2], result.parms_id(), result.scale(), pt);
    pms.eva->add_plain_inplace(result, pt);

    // cout << log2(abs(result.scale() - scale)) << endl;

    // 두 번째 항부터 반복
    for(int i=index-4; i>=0; i-=2)
    {
        // fi * (x^2) and mod reduce
        pms.encoder->encode(1.0, x_square.parms_id(), scale * correction_factors[cf_index++], pt);
        pms.eva->multiply_plain(x_square, pt, temp);
        pms.eva->rescale_to_next_inplace(temp);
        pms.eva->mod_reduce_to_inplace(temp, result.parms_id());

        // 곱하기
        pms.eva->multiply(result, temp, result);
        pms.eva->relinearize_inplace(result, pms.rlk);
        pms.eva->rescale_to_next_inplace(result);
        
        // 더하기
        pms.encoder->encode(coeff[i], result.parms_id(), result.scale(), pt);
        pms.eva->add_plain_inplace(result, pt);

        //debug
        // cout << log2(abs(result.scale() - scale)) << endl;
    }

    result.scale() = scale;
    return result;
}

// correction factor 계산. 짝수차 항만 존재하는 다항식에 해당.
vector<double> calculate_correction_factor(CKKS_params& pms, vector<double> coeff, Ciphertext& x)
{
    const auto& coeff_modulus = pms.context->get_context_data(x.parms_id())->parms().coeff_modulus();
    int coeff_modulus_index = coeff_modulus.size() - 1;
    int index = coeff.size() - 1;
    vector<double> correction_factors;

    double scale = x.scale();
    double scale_prime = (scale * scale) / static_cast<double_t>(coeff_modulus[coeff_modulus_index].value());

    Ciphertext temp = x;
    double scale_2 = static_cast<double_t>(coeff_modulus[coeff_modulus_index-1].value());
    double scale_i;

    correction_factors.push_back(scale_2 / scale_prime); // f1
    for(int i=1; i<index/2; i++)
    {
        scale_i = static_cast<double_t>(coeff_modulus[coeff_modulus_index-(i+1)].value());
        correction_factors.push_back((scale_2 * scale_i) / (scale * scale_prime));
        // correction_factors.push_back((scale_2) / (scale * scale_prime));
    }
    return correction_factors;
}

double calculate_cleanse_correction_factor(CKKS_params& pms, Ciphertext& x)
{
    const auto& coeff_modulus = pms.context->get_context_data(x.parms_id())->parms().coeff_modulus();
    int coeff_modulus_index = coeff_modulus.size() - 1;
    double correction_factor;
    
    double scale = x.scale();
    double scale_1 = static_cast<double_t>(coeff_modulus[coeff_modulus_index].value());
    double scale_2 = static_cast<double_t>(coeff_modulus[coeff_modulus_index-1].value());

    return (scale_1 * scale_1 * scale_2) / (scale * scale * scale);
}

Ciphertext cleanse(CKKS_params& pms, Ciphertext& x)
{
    double_t scale = x.scale();
    double corr_factor = calculate_cleanse_correction_factor(pms, x);

    // -2x^3 + 3x^2 = (-2x+3)x^2
    Ciphertext res, sqx;
    Plaintext pt;

    // (-2x+3) * f
    pms.encoder->encode(-2.0, x.parms_id(), 1.0, pt);
    pms.eva->multiply_plain(x, pt, res);
    pms.encoder->encode(3.0, res.parms_id(), res.scale(), pt);
    pms.eva->add_plain_inplace(res, pt);

    pms.encoder->encode(1.0, res.parms_id(), scale * corr_factor, pt);
    pms.eva->multiply_plain_inplace(res, pt);
    pms.eva->rescale_to_next_inplace(res);

    // x^2
    pms.eva->square(x, sqx);
    pms.eva->relinearize_inplace(sqx, pms.rlk);
    pms.eva->rescale_to_next_inplace(sqx);

    // (-2x+3)x^2
    pms.eva->multiply(res, sqx, res);
    pms.eva->relinearize_inplace(res, pms.rlk);
    pms.eva->rescale_to_next_inplace(res);


    // pms.encoder->encode(1.0, res.parms_id(), scale * corr_factor, pt);
    // pms.eva->multiply_plain_inplace(res, pt);

    // cout << log2(abs(res.scale() - scale)) << endl;
    res.scale() = scale;

    return res;
}

int evaluate_depth(vector<double> coeff, bool include_plain)
{
    int depth = coeff.size()/2;
    if(include_plain && coeff[coeff.size()-1] != 1.0)
        ++depth;
    return depth;
}

std::vector<double> linspace(double a, double b, std::size_t n) {
    std::vector<double> x;
    if (n == 0) return x;
    x.resize(n);

    if (n == 1) {
        x[0] = a; 
        return x;
    }

    double step = (b - a) / static_cast<double>(n - 1);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = a + step * static_cast<double>(i);
    }
    return x;
}


void compare_result(vector<double> res1, vector<double> res2, string title1, string title2, int size)
{
    cout << fixed << setprecision(6);
    vector<int> index = {0, 1, size-1};
    for(int i: index)
    {
        double log_diff = i == 0 ? log2(abs(1 - res2[i])) : log2(abs(res2[i]));
        cout << title1 << ": " << res1[i] << "\t\t" << title2 << ": " << res2[i] 
                << "\t\tlog2(diff): " << log_diff << endl;
    }
}

double compare_result_log(const std::vector<double>& res, int size, bool print)
{
    double min_val = 0;
    double max_val = -99;

    for (int i = 0; i < size; ++i) {
        double log_diff = i == 0 ? log2(abs(1 - res[i])) : log2(abs(res[i]));

        min_val = min(min_val, log_diff);
        max_val = max(max_val, log_diff);
    }
    if(print)
        std::cout << std::fixed << std::setprecision(6)
                << "log2(diff)"
                << " min: " << std::setw(9) << std::setfill('0') << min_val
                << ", max: " << std::setw(9) << std::setfill('0') << max_val
                << "\n";
    
    return max_val;
}

std::vector<double> exp_vector(const std::vector<double>& x, int exp)
{
    std::vector<double> res(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        res[i] = pow(x[i], exp);
    }
    return res;
}

std::vector<double> mult_vector(const std::vector<double>& a, const std::vector<double>& b)
{
    std::vector<double> res(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        res[i] = a[i] * b[i];
    }
    return res;
}

std::vector<double> add_vector(const std::vector<double>& a, const std::vector<double>& b)
{
    std::vector<double> res(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        res[i] = a[i] + b[i];
    }
    return res;
}

std::vector<double> mult_vector_plain(const std::vector<double>& a, double scalar)
{
    std::vector<double> res(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        res[i] = a[i] * scalar;
    }
    return res;
}

std::vector<double> add_vector_plain(const std::vector<double>& a, double scalar)
{
    std::vector<double> res(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        res[i] = a[i] + scalar;
    }
    return res;
}

std::vector<double> evaluate_function(const std::vector<double>& coeff, const std::vector<double>& x)
{
    std::vector<double> res(x.size(), 0.0);

    if (coeff.empty()) return res;

    for (std::size_t i = 0; i < x.size(); ++i) {
        double acc = 0.0;
        for (std::size_t j = coeff.size(); j-- > 0; ) {
            acc = acc * x[i] + coeff[j];
        }
        res[i] = acc;
    }
    return res;
}