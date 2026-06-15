#include "../include/evaluate.h"

void evaluate_helut(CKKS_params& pms, vector<double>& x, bool test_time)
{
    int x_bound = int(pms.x_bound);

    vector<double> real_result = x;
    vector<double> ckks_result;
    seal::Ciphertext ct = pms.encrypt(x);
    Plaintext pt;

    // 1. SqMethod - (1-2x^2/p^2)^(2^r)
    {
        // Base 연산 - +2 Depth
        pms.square_inplace(ct);
        pt = pms.encode((-2.0/(x_bound*x_bound)), ct);
        pms.mult_pt_ct_inplace(pt, ct, true);
        pt = pms.encode(1.0, ct);
        pms.add_pt_ct_inplace(pt, ct);

        //debug
        if(!test_time)
        {
            std::cout << "Base SqMethod." << std::endl;
            real_result = pt::mult_plain(real_result, (-2.0)/(x_bound*x_bound));
            real_result = pt::add_plain(real_result, 1);
            ckks_result = pms.decode_ctxt(ct);
            // compare_result_pt(ckks_result, real_result, x.size());
            // compare_result_IF(ckks_result, x.size());
            std::cout << "Remaining Levels: " << ct.coeff_modulus_size() << std::endl;
            std::cout << "---------------------------------" << std::endl;
        }

        // r번 제곱
        for(int i=0; i<pms.iter_HELUT.first; i++)
        {
            pms.square_inplace(ct);
            
            //debug
            if(!test_time)
            {
                std::cout << i+1 << "th SqMethod." << std::endl;
                real_result = pt::exp(real_result, 2);
                ckks_result = pms.decode_ctxt(ct);
                // compare_result_pt(ckks_result, real_result, x.size());
                // compare_result_IF(ckks_result, x.size());
                std::cout << "Remaining Levels: " << ct.coeff_modulus_size() << std::endl;
                std::cout << "---------------------------------" << std::endl;
            }
        }
    }


    // 2. Cleanse - -2x^3 + 3x^2
    vector<double> coeff_cleanse = {0.0, 0.0, 3.0, -2.0};
    seal::Ciphertext ct2;
    for(int i=0; i<pms.iter_HELUT.second; i++)
    {
        // x^2
        ct2 = pms.square(ct);
        
        // -2x+3
        pt = pms.encode(-2, ct, 1.0);
        pms.mult_pt_ct_inplace(pt, ct, false);
        pt = pms.encode(3, ct);
        pms.add_pt_ct_inplace(pt, ct);
        
        // x^2(-2x+3)
        // pms.eva->mod_reduce_to_next_inplace(ct);
        pt = pms.encode(1.0, ct);
        pms.mult_pt_ct_inplace(pt, ct, true);
        pms.mult_ct_ct_inplace(ct, ct2);

        if(!test_time)
        {
            real_result = pt::evaluate(coeff_cleanse, real_result);
            ckks_result = pms.decode_ctxt(ct);
            std::cout << "Cleanse " << i+1 << std::endl;
            // compare_result_pt(ckks_result, real_result, x.size());
            // compare_result_IF(ckks_result, x.size());
            std::cout << "Remaining Levels: " << ct.coeff_modulus_size() << std::endl;
            std::cout << "---------------------------------" << std::endl;
        }
    }

    // if(res_print)
    // {
    //     std::cout << "Final Result: " << std::endl;
    //     ckks_result = pms.decode_ctxt(ct);
    //     compare_result_pt(ckks_result, real_result, x.size());
    //     compare_result_IF(ckks_result, x.size());
    // }
}

void evaluate_multi_remez(CKKS_params& pms, vector<double>& x, std::vector<std::shared_ptr<Decomp>> dcmps, bool test_time)
{
    std::vector<double> real_result, ckks_result, coeff, corr_factors;
    real_result = x;
    seal::Ciphertext ct = pms.encrypt(x);

    // 3. 분해식 기반 함수 평가
    for(std::shared_ptr<Decomp> dcmp : dcmps)
    {        
        //debug
        if(!test_time)
        {
            std::cout << "Polynomial\n\t" << vec_to_str(dcmp->coeff) << std::endl;
            std::cout << "Decomp\n\t" << dcmp->restore_dcmp() << std::endl;
        }
        int start_level = ct.coeff_modulus_size();
        int target_use = ceil(log2(dcmp->coeff.size()));
        
        // routes, made_powers에 따라 x^i들을 구성
        std::map<int, seal::Ciphertext> powers = {{1, ct}};
        for(std::pair<int, int> route: dcmp->total_routes)
        {
            int a = route.first;
            int b = route.second;
            if(!dcmp->made_powers.contains(a+b))
                continue;
            seal::Ciphertext cta = powers.at(a);
            seal::Ciphertext ctb = powers.at(b);
            seal::Ciphertext ctres;
            pms.mult_ct_ct(cta, ctb, ctres);
            powers.insert({a+b, ctres});
        }
        
        // (recursive) x^i, p(x), q(x)
        double_t correction_factor = correction_factor_REMEZ(pms, dcmp->coeff, ct);
        seal::Ciphertext res = evaluate_polynomial(pms, ct, dcmp, powers, correction_factor);
        ct = res;
        
        int end_level = ct.coeff_modulus_size();
        if(target_use != (start_level - end_level))
        {
            std::cout << "LEVEL ERROR!, " << start_level << " -> " << end_level << std::endl;
            std::cout << "Target USE: " << target_use << std::endl;
            return;
        }
        //print
        if(!test_time)
        {
            real_result = pt::evaluate(dcmp->coeff, real_result);
            ckks_result = pms.decode_ctxt(ct);
            // compare_result_pt(ckks_result, real_result, x.size());
            // compare_result_IF(ckks_result, x.size());
            std::cout << std::endl;
            std::cout << "Remaining Levels: " << ct.coeff_modulus_size() << std::endl;
            std::cout << "---------------------------------" << std::endl;
        }
    }
}

double_t correction_factor_HELUT(CKKS_params &pms, seal::Ciphertext& x)
{
    const auto& coeff_modulus = pms.context->get_context_data(x.parms_id())->parms().coeff_modulus();
    int coeff_modulus_index = coeff_modulus.size() - 1;
    int iter_SqMethod = pms.iter_HELUT.first;
    int iter_Cleanse = pms.iter_HELUT.second;
    
    double_t correction_factor = 1.0;
    double_t scale = pms.scale;

    // CF for SqMethod.
    for(size_t i=0; i<1 + iter_SqMethod; i++)
    {
        correction_factor *= static_cast<double_t>(coeff_modulus[coeff_modulus_index].value()) / scale;
        --coeff_modulus_index;
    }

    //CF for Cleanse
    for(size_t i=0; i<2 * iter_Cleanse; i++)
    {
        correction_factor *= static_cast<double_t>(coeff_modulus[coeff_modulus_index].value()) / scale;
        --coeff_modulus_index;
    }
    return correction_factor;
}

double_t correction_factor_REMEZ(CKKS_params& pms, std::vector<double>& coeff, seal::Ciphertext& ct)
{
    int required_depth = ceil(log2((int)coeff.size()));
    const auto& coeff_modulus = pms.context->get_context_data(ct.parms_id())->parms().coeff_modulus();
    
    double_t correction_factor = 1.0;
    double_t target_scale = ct.scale();
    int mod_index = ct.coeff_modulus_size() - 1;
    for(int i=0; i<required_depth; ++i)
    {
        correction_factor *= (static_cast<double_t>(coeff_modulus[mod_index--].value()) / target_scale);
    }
    return correction_factor;
}

seal::Ciphertext evaluate_xi(CKKS_params& pms, seal::Ciphertext& ct, XI xi, std::map<int, seal::Ciphertext> powers, double_t correction_factor)
{
    // 1. x^i
    if(!xi.multA)
    {
        return powers.at(xi.n);
    }
    else
    {
        seal::Ciphertext ct_xi;
        seal::Plaintext pt;
        
        if(!powers.contains(xi.n))
        {
            seal::Ciphertext cta, ctb;
            std::pair<int, int> final_route = xi.route.back();
            cta = powers.at(final_route.first);
            pt = pms.encode(xi.coeff * correction_factor, cta);  
            pms.mult_pt_ct_inplace(pt, cta, true);
            ctb = powers.at(final_route.second);
            pms.mult_ct_ct(cta, ctb, ct_xi);
        }
        else
        {
            ct_xi = powers.at(xi.n);
            pt = pms.encode(xi.coeff * correction_factor, ct_xi);        
            pms.mult_pt_ct_inplace(pt, ct_xi, true);
        }
        return ct_xi;
    }
}

seal::Ciphertext evaluate_polynomial(CKKS_params& pms, seal::Ciphertext& ct, std::shared_ptr<Decomp> dcmp, std::map<int, seal::Ciphertext> powers, double_t correction_factor)
{
    seal::Ciphertext res;
    bool flag_p, flag_q;
    flag_p = dcmp->dcmp_p ? true : false;
    flag_q = dcmp->dcmp_q ? true : false;

    // Leaf node
    if(!flag_p && !flag_q)
    {
        bool flag = false;
        seal::Ciphertext x;
        seal::Plaintext pt;
        // 상수항
        //ax, a+bx^2, a+x^2
        if(dcmp->coeff[0] != 0)
        {
            seal::Ciphertext ct2 = pms.encrypt(dcmp->coeff[0], ct);
            if(!flag)
            {
                res = ct2;
                flag = true;
            }
        }
        for(int i=1; i<dcmp->coeff.size(); i++)
        {
            if(dcmp->coeff[i] == 0.0)
                continue;
            x = powers.at(i);
            double scale = dcmp->coeff[i] == 1 ? 1.0 : x.scale();
            bool rescale = dcmp->coeff[i] == 1 ? false : true;
            pt = pms.encode(dcmp->coeff[i], x, scale);
            pms.mult_pt_ct_inplace(pt, x, rescale);
            if(!flag)
            {
                res = x;
                flag = true;
            }
            else
            {
                pms.add_ct_ct(res, x, res);
            }
        }
        return res;
    }

    seal::Ciphertext res_xi, res_px, res_qx;
    if(flag_p)
    {
        res_px = evaluate_polynomial(pms, ct, dcmp->dcmp_p, powers);
        res_xi = evaluate_xi(pms, ct, dcmp->xi, powers, correction_factor);
        pms.mult_ct_ct(res_xi, res_px, res);
    }
    if(flag_q)
    {
        res_qx = evaluate_polynomial(pms, ct, dcmp->dcmp_q, powers);
        pms.add_ct_ct(res, res_qx, res);
    }     
    
    return res;
}

/*
    // // 1. coeff 1차
    // coeff = coeffs[0];
    // corr_factors = calculate_correction_factor(pms, coeff, ct);
    // ct = eval_even(pms, coeff, corr_factors, ct);
    // if(res_print)
    // {
    //     real_result = evaluate_function(coeff, x);
    //     if(mid_print)
    //     {
    //         ckks_result = pms.decode_ctxt(ct);
    //         cout << "Coeff Evaluation 0, degree " << coeff.size()-1 << endl;
    //         // compare_result(real_result, ckks_result, "REAL", "CKKS", x.size());
    //         double max_err = compare_result_log(ckks_result, x.size(), true);
    //         cout << "Remaining Levels: " << ct.coeff_modulus_size() << endl;
    //         cout << "------------\n";
    //     }
    // }
    // x = real_result;

    // // 2. coeff 2차~
    // for(int i=1; i<mi.n - mi.s; i++)
    // {
    //     coeff = coeffs[i];
    //     corr_factor = calculate_correction_factor_odd(pms, coeff, ct);
    //     ct = eval_odd(pms, coeff, corr_factor, ct);
    //     if(res_print)
    //     {
    //         real_result = evaluate_function(coeff, x);
    //         if(mid_print)
    //         {
    //             ckks_result = pms.decode_ctxt(ct);
    //             cout << "Coeff Evaluation " << i+1 << ", degree " << coeff.size()-1 << endl;
    //             // compare_result(real_result, ckks_result, "REAL", "CKKS", x.size());
    //             double max_err = compare_result_log(ckks_result, x.size(), true);
    //             cout << "Remaining Levels: " << ct.coeff_modulus_size() << endl;
    //             cout << "------------\n";
    //         }
    //     }
    //     x = real_result;
    // }

    // // 3. Cleanse
    // vector<double> coeff_cleanse = {0.0, 0.0, 3.0, -2.0};
    // for(int i=0; i<mi.s; i++)
    // {
    //     ct = cleanse(pms, ct);
    //     if(res_print)
    //     {
    //         real_result = evaluate_function(coeff_cleanse, real_result);
    //         if(mid_print)
    //         {
    //             ckks_result = pms.decode_ctxt(ct);
    //             cout << "Cleanse result" << endl;
    //             // compare_result(real_result, ckks_result, "REAL", "CKKS", x.size());
    //             double max_err = compare_result_log(ckks_result, x.size(), true);
    //             cout << "------------\n";
    //         }
    //     }
    // }

    // if(res_print)
    // {
    //     cout << "Final Result: " << endl;
    //     ckks_result = pms.decode_ctxt(ct);
    //     compare_result(real_result, ckks_result, "REAL", "CKKS", real_result.size());
    //     // compare_result_log(ckks_result, x.size(), true);
    // }

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

Ciphertext eval_even(CKKS_params& pms, vector<double> coeff, vector<double> correction_factors, Ciphertext& x)
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

Ciphertext eval_odd(CKKS_params& pms, vector<double> coeff, double correction_factor, Ciphertext& x)
{
    int max_deg = coeff.size() - 1;
    int cf_index = 0;
    double scale = x.scale();

    // 1. 최고차항의 계수가 1이 되도록 계수 재정의
    double outer_coeff = coeff[max_deg];
    for(int i=1; i<=max_deg; i+=2)
        coeff[i] /= outer_coeff;

    // 2. 항 계산
    //식 형태는 (짝수차 함수) * ax
    Plaintext pt;
    Ciphertext x_square, temp, result;

    //x^2
    pms.eva->square(x, x_square);
    pms.eva->relinearize_inplace(x_square, pms.rlk);
    pms.eva->rescale_to_next_inplace(x_square);

    // 첫 번째 항 = x^2 + b/a
    pms.encoder->encode(coeff[max_deg-2], x_square.parms_id(), x_square.scale(), pt);
    pms.eva->add_plain(x_square, pt, result);

    // cout << log2(abs(result.scale() - scale)) << endl;

    // 두 번째 항
    if(max_deg == 5)
    {
        // 곱하기
        pms.eva->multiply(result, x_square, result);
        pms.eva->relinearize_inplace(result, pms.rlk);
        pms.eva->rescale_to_next_inplace(result);
        
        // 더하기
        pms.encoder->encode(coeff[1], result.parms_id(), result.scale(), pt);
        pms.eva->add_plain_inplace(result, pt);
    }

    //마지막 항(ax)
    pms.encoder->encode(outer_coeff, x.parms_id(), scale * correction_factor, pt);
    pms.eva->multiply_plain(x, pt, temp);
    pms.eva->rescale_to_next_inplace(temp);
    pms.eva->mod_reduce_to_inplace(temp, result.parms_id());

    pms.eva->multiply(result, temp, result);
    pms.eva->relinearize_inplace(result, pms.rlk);
    pms.eva->rescale_to_next_inplace(result);

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

// correction factor 계산. 홀수차 다항식용.
double calculate_correction_factor_odd(CKKS_params& pms, vector<double> coeff, Ciphertext& x)
{
    const auto& coeff_modulus = pms.context->get_context_data(x.parms_id())->parms().coeff_modulus();
    int coeff_modulus_index = coeff_modulus.size() - 1;
    int max_deg = coeff.size() - 1;

    double scale, scale_p, scale_pp, scale_ppp;
    double f;

    scale = x.scale();
    scale_p = static_cast<double_t>(coeff_modulus[coeff_modulus_index].value());
    scale_pp = static_cast<double_t>(coeff_modulus[coeff_modulus_index-1].value());

    if(max_deg == 3)
        f = (scale_p * scale_p * scale_pp) / (scale * scale * scale);
    else if(max_deg == 5)
    {
        scale_ppp = static_cast<double_t>(coeff_modulus[coeff_modulus_index-2].value());
        f = (pow(scale_p, 3) * scale_pp * scale_ppp) / pow(scale, 5);
    }
    return f;
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

// 다항식 타입 체크
int check_coeff_type(vector<double> coeff)
{
    int max_deg = coeff.size() - 1;
    if(max_deg % 2 == 0)
    {
        for(int i=1; i<=max_deg; i+=2)
            if(coeff[i] != 0.0)
                return 2;
        return 0;
    }
    else if(max_deg % 2 == 1)
    {
        for(int i=0; i<=max_deg; i+=2)
            if(coeff[i] != 0.0)
                return 2;
        return 1;
    }
    return 2;
}

int evaluate_depth(vector<double> coeff, bool include_plain)
{
    // 1. 다항식 타입(all, odd, even) 구별
    int coeff_type = check_coeff_type(coeff);

    // 2. 다항식 타입에 따른 depth반환
    switch(coeff_type)
    {
        case 0: // 짝함수
            return coeff.size() / 2 + 1;
        case 1: // 홀함수
            return coeff.size() / 2;
        case 2: // 일반함수
            break;
    }
    return -99;
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
*/

// 분해식 정보 기반 암호문 평가
/* CT eval_poly_ct(CKKS_params& pms, shared_ptr<Decomp> dcmp, CT x, ErrBound& eb, int copy_count, bool step_debug, bool res_debug)
{
    EvalStep es(dcmp);
    // es.print_step();

    //plain value
    std::map<std::string, double> coeffs;
    std::vector<double> ptx = std::get<0>(x);
    std::map<std::string, std::vector<double>> ppowers;
    std::map<std::string, std::vector<double>> pterms;
    ppowers.insert_or_assign("P1", ptx);
    std::vector<double> pvalue1, pvalue2, padd_res, pmult_res;

    //seal::Ciphertext datas
    seal::Ciphertext ctx = std::get<1>(x);
    std::map<std::string, seal::Ciphertext> cpowers;
    std::map<std::string, seal::Ciphertext> cterms;
    cpowers.insert_or_assign("P1", ctx); // Assert P1
    seal::Ciphertext cvalue1, cvalue2, cadd_res, cmult_res;
    seal::Plaintext pt;

    //fct::Ciphertext datas
    fct::Ciphertext fctx = std::get<2>(x);
    std::map<std::string, fct::Ciphertext> fpowers;
    std::map<std::string, fct::Ciphertext> fterms;
    fpowers.insert_or_assign("P1", fctx);
    fct::Ciphertext fvalue1, fvalue2, fadd_res, fmult_res;
    
    // // Correction factor 연산
    // double correction_factor = 1.0;
    // int required_depth = ceil(log2((int)dcmp->coeff.size()));
    // const auto& coeff_modulus = pms.context->get_context_data(ctx.parms_id())->parms().coeff_modulus();

    // int mod_index = ctx.coeff_modulus_size() - 1;
    // double previous_scale = ctx.scale();
    // double next_scale;
    // for(int i=0; i<required_depth; ++i)
    // {
    //     next_scale = (previous_scale * previous_scale) / coeff_modulus[mod_index--].value();
    //     correction_factor = next_scale / previous_scale;
    // }

    // int apply_index = 0;
    // for(int i=0; i<es.eval_step.size()-1; i++)
    // {
    //     auto step = es.eval_step[i];
    //     auto nstep = es.eval_step[i+1];
    //     if(step.op == 'o' && nstep.op == 'x')
    //         apply_index = i;
    // }

    // 평가
    int step_count = 0;
    std::tuple<std::vector<double>, seal::Ciphertext, fct::Ciphertext> step_res;
    for(auto step: es.eval_step)
    {
        if(step_debug)
        {
            std::cout << "Operation" << std::endl;
            es.print_line(step, step_count);
        }
        switch(step.op)
        {
            case 'o':
                if(step.key1[0] != 'C')
                {
                    if(step.key1[0] == 'P')
                        continue;
                    throw std::invalid_argument(std::format("Main::RES: invalid key data {} {}", step.key1, step.key2));
                }
                coeffs.insert_or_assign(step.key1, std::stod(step.key2));
                break;
            case '+':
                if(step.key1[0] == 'T' && step.key2[0] == 'C')
                {
                    cvalue1 = cterms.at(step.key1);
                    pt = pms.encode(coeffs.at(step.key2), cvalue1);
                    pms.add_pt_ct(pt, cvalue1, cadd_res);

                    padd_res = pt::add_plain(pterms.at(step.key1), coeffs.at(step.key2));
                    fadd_res = fct::add_pt(fterms.at(step.key1), coeffs.at(step.key2));
                }
                else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                {
                    cvalue1 = cterms.at(step.key1);
                    cvalue2 = cterms.at(step.key2);
                    pms.add_ct_ct(cvalue1, cvalue2, cadd_res);

                    padd_res = pt::add(pterms.at(step.key1), pterms.at(step.key2));
                    fadd_res = fct::add(fterms.at(step.key1), fterms.at(step.key2));
                }
                else if(step.key1[0] == 'P' && step.key2[0] == 'C')
                {
                    cvalue1 = cpowers.at(step.key1);
                    pt = pms.encode(coeffs.at(step.key2), cvalue1);
                    pms.add_pt_ct(pt, cvalue1, cadd_res);

                    padd_res = pt::add_plain(ppowers.at(step.key1), coeffs.at(step.key2));
                    fadd_res = fct::add_pt(fpowers.at(step.key1), coeffs.at(step.key2));
                }
                else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                {
                    cvalue1 = cpowers.at(step.key1);
                    cvalue2 = cterms.at(step.key2);
                    pms.add_ct_ct(cvalue1, cvalue2, cadd_res);

                    padd_res = pt::add(ppowers.at(step.key1), pterms.at(step.key2));
                    fadd_res = fct::add(fpowers.at(step.key1), fterms.at(step.key2));
                }
                else
                {
                    throw std::invalid_argument(std::format("Main::ADD: invalid key data {} {}", step.key1, step.key2));
                }
                pterms.insert_or_assign(step.save_key, padd_res);
                cterms.insert_or_assign(step.save_key, cadd_res);
                fterms.insert_or_assign(step.save_key, fadd_res);
                step_res = {padd_res, cadd_res, fadd_res};
                break;

            case 'x':
                if(step.key1[0] == 'P' && step.key2[0] == 'P')
                {
                    cvalue1 = cpowers.at(step.key1);
                    cvalue2 = cpowers.at(step.key2);
                    pms.mult_ct_ct(cvalue1, cvalue2, cmult_res);
                    

                    pmult_res = pt::mult(ppowers.at(step.key1), ppowers.at(step.key2));
                    if(step.key1 == step.key2)
                        fmult_res = fct::square(fpowers.at(step.key1), eb);
                    else
                        fmult_res = fct::mult(fpowers.at(step.key1), fpowers.at(step.key2), eb);
                }
                else if(step.key1[0] == 'C' && step.key2[0] == 'P')
                {
                    double coeff = coeffs.at(step.key1);
                    cvalue1 = cpowers.at(step.key2);
                    
                    // if(apply_index == step_count)
                    //     coeff *= correction_factor;
                    pt = pms.encode(coeff, cvalue1);
                    pms.mult_pt_ct(pt, cvalue1, cmult_res, true);

                    pmult_res = pt::mult_plain(ppowers.at(step.key2), coeffs.at(step.key1));
                    fmult_res = fct::mult_pt(coeffs.at(step.key1), fpowers.at(step.key2), eb);
                }
                else if(step.key1[0] == 'P' && step.key2[0] == 'C')
                {
                    cvalue1 = cpowers.at(step.key1);
                    double coeff = coeffs.at(step.key2);
                    // if(apply_index == step_count)
                    //     coeff *= correction_factor;
                    pt = pms.encode(coeff, cvalue1);
                    pms.mult_pt_ct(pt, cvalue1, cmult_res, true);

                    pmult_res = pt::mult_plain(ppowers.at(step.key1), coeffs.at(step.key2));
                    fmult_res = fct::mult_pt(coeffs.at(step.key2), fpowers.at(step.key1), eb);
                }
                else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                {
                    cvalue1 = cpowers.at(step.key1);
                    cvalue2 = cterms.at(step.key2);
                    pms.mult_ct_ct(cvalue1, cvalue2, cmult_res);

                    pmult_res = pt::mult(ppowers.at(step.key1), pterms.at(step.key2));
                    fmult_res = fct::mult(fpowers.at(step.key1), fterms.at(step.key2), eb);
                }
                else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                {
                    cvalue1 = cterms.at(step.key1);
                    cvalue2 = cterms.at(step.key2);
                    pms.mult_ct_ct(cvalue1, cvalue2, cmult_res);

                    pmult_res = pt::mult(pterms.at(step.key1), pterms.at(step.key2));
                    fmult_res = fct::mult(fterms.at(step.key1), fterms.at(step.key2), eb);
                }
                else if(step.key1[0] == 'T' && step.key2[0] == 'P')
                {
                    cvalue1 = cterms.at(step.key1);
                    cvalue2 = cpowers.at(step.key2);
                    pms.mult_ct_ct(cvalue1, cvalue2, cmult_res);

                    pmult_res = pt::mult(pterms.at(step.key1), ppowers.at(step.key2));
                    fmult_res = fct::mult(fterms.at(step.key1), fpowers.at(step.key2), eb);
                }
                else
                {
                    throw std::invalid_argument(std::format("Main::MUL: invalid key data {} {}", step.key1, step.key2));
                }
                switch(step.save_key[0])
                {
                    case 'P':
                        ppowers.insert_or_assign(step.save_key, pmult_res);
                        cpowers.insert_or_assign(step.save_key, cmult_res);
                        fpowers.insert_or_assign(step.save_key, fmult_res);
                        break;
                    case 'T':
                        pterms.insert_or_assign(step.save_key, pmult_res);
                        cterms.insert_or_assign(step.save_key, cmult_res);
                        fterms.insert_or_assign(step.save_key, fmult_res);
                        break;
                    default:
                        break;
                }
                step_res = {pmult_res, cmult_res, fmult_res};
                break;
            default:
                break;
        }
        step_count++;
        
        // debug - compare at each step    
        if(step_debug && step.op != 'o')    
        {
            compare_boundary(step_res, pms, copy_count);
            std::cout << "-----------------------------------------------------" << std::endl;
        }
    }
    // std::vector<double> pres = pterms.at(es.eval_step.back().save_key);
    // seal::Ciphertext cres = cterms.at(es.eval_step.back().save_key);
    // fct::Ciphertext fres = fterms.at(es.eval_step.back().save_key);
    // CT res = {pres, cres, fres};
    if(res_debug)
    {
        compare_boundary(step_res, pms, copy_count);
    }
    return step_res;
} */

// 분해식 정보 기반 seal::Ciphertext 평가

seal::Ciphertext ct::evaluate(seal::Ciphertext& x, EvalStep es, CKKS_params& pms)
{
    //seal::Ciphertext datas
    seal::Ciphertext ctx = x;
    std::map<std::string, double> coeffs;
    std::map<std::string, seal::Ciphertext> cpowers;
    std::map<std::string, seal::Ciphertext> cterms;
    cpowers.insert_or_assign("P1", ctx); // Assert P1
    seal::Ciphertext cvalue1, cvalue2, cadd_res, cmult_res;
    seal::Plaintext pt;

    // 평가
    int step_count = 0;
    seal::Ciphertext step_res;
    for(auto step: es.eval_step)
    {
        switch(step.op)
        {
            case 'o':
                if(step.key1[0] != 'C')
                {
                    if(step.key1[0] == 'P')
                        continue;
                    throw std::invalid_argument(std::format("Main::RES: invalid key data {} {}", step.key1, step.key2));
                }
                coeffs.insert_or_assign(step.key1, std::stod(step.key2));
                break;
            case '+':
                if(step.key1[0] == 'T' && step.key2[0] == 'C')
                {
                    cvalue1 = cterms.at(step.key1);
                    pt = pms.encode(coeffs.at(step.key2), cvalue1);
                    pms.add_pt_ct(pt, cvalue1, cadd_res);
                }
                else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                {
                    cvalue1 = cterms.at(step.key1);
                    cvalue2 = cterms.at(step.key2);
                    pms.add_ct_ct(cvalue1, cvalue2, cadd_res);
                }
                else if(step.key1[0] == 'P' && step.key2[0] == 'C')
                {
                    cvalue1 = cpowers.at(step.key1);
                    pt = pms.encode(coeffs.at(step.key2), cvalue1);
                    pms.add_pt_ct(pt, cvalue1, cadd_res);
                }
                else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                {
                    cvalue1 = cpowers.at(step.key1);
                    cvalue2 = cterms.at(step.key2);
                    pms.add_ct_ct(cvalue1, cvalue2, cadd_res);
                }
                else
                {
                    throw std::invalid_argument(std::format("Main::ADD: invalid key data {} {}", step.key1, step.key2));
                }
                cterms.insert_or_assign(step.save_key, cadd_res);
                step_res = cadd_res;
                break;

            case 'x':
                if(step.key1[0] == 'P' && step.key2[0] == 'P')
                {
                    cvalue1 = cpowers.at(step.key1);
                    cvalue2 = cpowers.at(step.key2);
                    pms.mult_ct_ct(cvalue1, cvalue2, cmult_res);
                }
                else if(step.key1[0] == 'C' && step.key2[0] == 'P')
                {
                    double coeff = coeffs.at(step.key1);
                    cvalue1 = cpowers.at(step.key2);
                    pt = pms.encode(coeff, cvalue1);
                    pms.mult_pt_ct(pt, cvalue1, cmult_res, true);
                }
                else if(step.key1[0] == 'P' && step.key2[0] == 'C')
                {
                    cvalue1 = cpowers.at(step.key1);
                    double coeff = coeffs.at(step.key2);
                    pt = pms.encode(coeff, cvalue1);
                    pms.mult_pt_ct(pt, cvalue1, cmult_res, true);
                }
                else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                {
                    cvalue1 = cpowers.at(step.key1);
                    cvalue2 = cterms.at(step.key2);
                    pms.mult_ct_ct(cvalue1, cvalue2, cmult_res);
                }
                else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                {
                    cvalue1 = cterms.at(step.key1);
                    cvalue2 = cterms.at(step.key2);
                    pms.mult_ct_ct(cvalue1, cvalue2, cmult_res);
                }
                else if(step.key1[0] == 'T' && step.key2[0] == 'P')
                {
                    cvalue1 = cterms.at(step.key1);
                    cvalue2 = cpowers.at(step.key2);
                    pms.mult_ct_ct(cvalue1, cvalue2, cmult_res);
                }
                else
                {
                    throw std::invalid_argument(std::format("Main::MUL: invalid key data {} {}", step.key1, step.key2));
                }
                switch(step.save_key[0])
                {
                    case 'P':
                        cpowers.insert_or_assign(step.save_key, cmult_res);
                        break;
                    case 'T':
                        cterms.insert_or_assign(step.save_key, cmult_res);
                        break;
                    default:
                        break;
                }
                step_res = cmult_res;
                break;
            default:
                break;
        }
        step_count++;
    }
    return step_res;
}

// Cleanse함수 평가
/* CT eval_poly_cl(CKKS_params& pms, CT x, ErrBound& eb, int copy_count, bool step_debug, bool res_debug)
{
    //plain value
    std::map<std::string, double> coeffs;
    std::vector<double> ptx = std::get<0>(x);
    std::map<std::string, std::vector<double>> ppowers;
    std::map<std::string, std::vector<double>> pterms;
    ppowers.insert_or_assign("P1", ptx);
    std::vector<double> pvalue1, pvalue2, padd_res, pmult_res;

    //seal::Ciphertext datas
    seal::Ciphertext ctx = std::get<1>(x);
    std::map<std::string, seal::Ciphertext> cpowers;
    std::map<std::string, seal::Ciphertext> cterms;
    cpowers.insert_or_assign("P1", ctx); // Assert P1
    seal::Ciphertext cvalue1, cvalue2, cadd_res, cmult_res;
    seal::Plaintext pt;

    //fct::Ciphertext datas
    fct::Ciphertext fctx = std::get<2>(x);
    std::map<std::string, fct::Ciphertext> fpowers;
    std::map<std::string, fct::Ciphertext> fterms;
    fpowers.insert_or_assign("P1", fctx);
    fct::Ciphertext fvalue1, fvalue2, fadd_res, fmult_res;

    // x^2
    std::vector<double> ptx2 = pt::mult(ptx, ptx);
    seal::Ciphertext ctx2 = pms.square(ctx);
    fct::Ciphertext fctx2 = fct::square(fctx, eb);
    // debug - compare at each step
    if(step_debug)
    {
        compare_boundary({ptx2, ctx2, fctx2}, pms, copy_count);
        std::cout << "-----------------------------------------------------" << std::endl;
    }

    // -2x+3
    std::vector<double> ptx3 = pt::mult_plain(ptx2, -2);
    ptx3 = pt::add_plain(ptx3, 3);
    pt = pms.encode(-2, ctx2, 1.0);
    seal::Ciphertext ctx3;
    pms.mult_pt_ct(pt, ctx2, ctx3, false);
    pt = pms.encode(3, ctx3);
    pms.add_pt_ct_inplace(pt, ctx3);
    fct::Ciphertext fctx3 = fct::mult_pt(-2, fctx2, eb, false);
    fctx3 = fct::add_pt(fctx3, 3.0);
    // debug - compare at each step    
    if(step_debug)    
    {
        compare_boundary({ptx3, ctx3, fctx3}, pms, copy_count);
        std::cout << "-----------------------------------------------------" << std::endl;
    }

    // x^2(-2x+3)
    std::vector<double> pt_res = pt::mult(ptx2, ptx3);
    seal::Ciphertext ct_res;
    pms.mult_ct_ct(ctx2, ctx3, ct_res, false);
    fct::Ciphertext fct_res = fct::mult(fctx2, fctx3, eb);
    CT res = {pt_res, ct_res, fct_res};
    // debug - compare at each step    
    if(step_debug)    
    {
        compare_boundary(res, pms, copy_count);
        std::cout << "-----------------------------------------------------" << std::endl;
    }

    // debug - Final Result
    if(res_debug)    
    {
        compare_boundary(res, pms, copy_count);
        // std::cout << "-----------------------------------------------------" << std::endl;
    }
    return res;
}
 */


// 비교용
void compare_boundary(std::vector<double> ptx, seal::Ciphertext ctx, std::pair<double, double> fct_interval, CKKS_params& pms, int copy_count)
{            

    std::vector<double> pt_res = ptx;
    std::vector<double> ct_res = pms.decode_ctxt(ctx);
    double fct_low = fct_interval.first;
    double fct_high = fct_interval.second;

    bool boundary_flag;

    int answer = std::round(pt_res[0]);
    double ct_high = *std::max_element(ct_res.begin(), ct_res.end());
    double ct_low = *std::min_element(ct_res.begin(), ct_res.end());

    // Boundary check
    boundary_flag = fct_low <= ct_low && ct_high <= fct_high ? true : false;

    // print result
    if(!boundary_flag)
    {
        std::cout 
            << std::format("Boundary FAIL in x={}\n", 0)
            << std::format("\tAnswer\t{}\n", clearNum(answer))
            << std::format("\tHigh boundary\t{}\n", clearNum(fct_high))
            << std::format("\tDecrypt result high\t{}\n", clearNum(ct_high))
            << std::format("\tDecrypt result low\t{}\n", clearNum(ct_low))
            << std::format("\tLow boundary\t{}\n", clearNum(fct_low));
    }
    if(boundary_flag)
    {
        double margin = std::max(std::abs(fct_high - ct_high), std::abs(fct_low - ct_low));
        std::cout << std::format("Boundary PASS. margin: {}", clearNum(margin)) << std::endl;
    }
}

// Assume all data are target to equal value.
double compare_precision(std::vector<double> ptx, seal::Ciphertext ctx, std::pair<double, double> fct_interval, CKKS_params& pms, int copy_count)
{
    std::vector<double> pt_res = ptx;
    std::vector<double> ct_res = pms.decode_ctxt(ctx);
    double fct_low = fct_interval.first;
    double fct_high = fct_interval.second;

    // 최종 결과 비교 - 목표 정밀도에 부합하는가?
    double precision = pms.target_precision;
    bool boundary_flag, precision_flag;

    int answer = std::round(pt_res[0]);
    // double fct_high = *std::max_element(fct_res_high.begin(), fct_res_high.end());
    double ct_high = *std::max_element(ct_res.begin(), ct_res.end());
    double ct_low = *std::min_element(ct_res.begin(), ct_res.end());
    // double fct_low = *std::min_element(fct_res_low.begin(), fct_res_low.end());

    // Boundary check
    boundary_flag = fct_low <= ct_low && ct_high <= fct_high ? true : false;
    // Precision check
    double diff = std::max(std::abs(ct_high - answer), std::abs(answer - ct_low));
    precision_flag = diff <= precision ? true : false;

    // print result
    if(!boundary_flag && !precision_flag)
    {
        std::cout 
            << std::format("Boundary FAIL and Precision FAIL in x={}\n", 0)
            << std::format("\tAnswer\t{}\n", clearNum(answer))
            << std::format("\tHigh boundary\t{}\n", clearNum(fct_high))
            << std::format("\tDecrypt result high\t{}\n", clearNum(ct_high))
            << std::format("\tDecrypt result low\t{}\n", clearNum(ct_low))
            << std::format("\tLow boundary\t{}\n", clearNum(fct_low));
    }
    else if(!boundary_flag && precision_flag)
    {
        std::cout
            << std::format("Boundary FAIL, but Precision acheived in x={}\n", 0)
            << std::format("\tAnswer\t{}\n", clearNum(answer))
            << std::format("\tHigh boundary\t{}\n", clearNum(fct_high))
            << std::format("\tDecrypt result high\t{}\n", clearNum(ct_high))
            << std::format("\tDecrypt result low\t{}\n", clearNum(ct_low))
            << std::format("\tLow boundary\t{}\n", clearNum(fct_low));
    }
    else if(boundary_flag && !precision_flag)
    {
        std::cout
            << std::format("Boundary PASS, but Precision failed in x={}\n", 0)
            << std::format("\tAnswer\t{}\n", clearNum(answer))
            << std::format("\tHigh boundary\t{}\n", clearNum(fct_high))
            << std::format("\tDecrypt result high\t{}\n", clearNum(ct_high))
            << std::format("\tDecrypt result low\t{}\n", clearNum(ct_low))
            << std::format("\tLow boundary\t{}\n", clearNum(fct_low));
    }

    /* for(int i=0; i<ptx.size() / copy_count; i++)
    {
        boundary_flag = true;
        precision_flag = true;
        double pt_target = (ptx[i] == 0) ? 1.0 : 0.0;

        for(int j=0; j<copy_count; j++)
        {
            int idx = (i*copy_count) + j;
            double low = fct_res_low[idx];
            double high = fct_res_high[idx];
            double dec_value = ct_res[idx];
            // fct boundary check
            if(dec_value < low || high < dec_value)
            {
                boundary_flag = false;
                final_flag = false;
            }
            
            // target precision check
            double diff = std::abs(dec_value - pt_target);
            if(diff > precision)
            {
                precision_flag = false;
                final_flag = false;
                std::cout 
                    << std::format("Target Precision FAIL in x={}\n", i)
                    << std::format("\tDecrypt result\t{}\n", clearNum(dec_value));

                if(dec_value > (pt_target + precision))
                    std::cout << std::format("\tHigh boundary\t{}\n", clearNum(pt_target+precision));
                else if(dec_value < (pt_target - precision))
                    std::cout << std::format("\tLow boundary\t{}\n", clearNum(pt_target-precision));
            }

            if(!boundary_flag && !precision_flag)
            {
                std::cout 
                    << std::format("Boundary FAIL in x={}\n", i)
                    << std::format("\tLow Boundary\t{}\n", clearNum(low))
                    << std::format("\tDecrypt result\t{}\n", clearNum(dec_value))
                    << std::format("\tHigh boundary\t{}\n", clearNum(high));
                break;
            }
            else if(!boundary_flag && precision_flag)
            {
                std::cout
                    << std::format("Boundary FAIL, but Precision acheived in x={}\n", i)
                    << std::format("\tLow Boundary\t{}\n", clearNum(low))
                    << std::format("\tDecrypt result\t{}\n", clearNum(dec_value))
                    << std::format("\tHigh boundary\t{}\n", clearNum(high));
            }
            else if(boundary_flag && !precision_flag)
            {
                std::cout
                    << std::format("Boundary PASS, but Precision failed in x={}\n", i)
                    << std::format("\tLow Boundary\t{}\n", clearNum(low))
                    << std::format("\tDecrypt result\t{}\n", clearNum(dec_value))
                    << std::format("\tHigh boundary\t{}\n", clearNum(high));
            }
            else 
            {
                max_margin = std::max(max_margin, diff);
            }
        } */
    
    if(boundary_flag && precision_flag)
    {
        double margin = std::max(std::abs(fct_high - ct_high), std::abs(fct_low - ct_low));
        std::cout << std::format("All result PASS. margin: {}", clearNum(margin)) << std::endl;
        // std::cout << std::format("\t{} < {}", clearNum(ct_high), clearNum(fct_high)) << std::endl;
        // std::cout << std::format("\t{} < {}", clearNum(fct_low), clearNum(ct_low)) << std::endl;
        return margin;
    }
    return -99;
}

// Compare result of ciphertext with plaintext value.
/* void compare_result_pt(const std::vector<double>& ct_res, const std::vector<double>& pt_res, int size)
{
    std::cout << fixed << setprecision(6);
    double min_err = 99.0;
    double max_err = -99.0;
    int min_index, max_index;
    double log_diff;

    for(int i=0; i<size; i++)
    {
        log_diff = log2(ct_res[i] - pt_res[i]);
        if(log_diff > max_err)
        {
            max_err = log_diff;
            max_index = i;
        }
        else if(log_diff < min_err)
        {
            min_err = log_diff;
            min_index = i;
        }
    }
    std::cout << "Comparison with plaintext value" << std::endl;
    std::cout << "\tmax_err=" << max_err << std::endl;
    std::cout << "\tmin_err=" << min_err << std::endl;
} */

// Compare result of ciphertext with target value.
/* void compare_result_IF(const std::vector<double>& res, int size)
{
    double min_err = 99.0;
    double max_err = -99.0;
    double log_diff;
    for (int i = 0; i < size; ++i) {
        log_diff = log2(abs(round(res[i]) - res[i]));
        max_err = log_diff > max_err ? log_diff : max_err;
        min_err = log_diff < min_err ? log_diff : min_err;
    }
    std::cout << "Comparison with target value" << std::endl;
    std::cout << "\tmax_err=" << max_err << std::endl;
    std::cout << "\tmin_err=" << min_err << std::endl;
} */


/* std::vector<double> evaluate_polynomial_dcmp(const std::vector<double>& coeff, const std::vector<double>& x, EvalStep es)
{
    //plain value
    std::map<std::string, double> coeffs;
    std::map<std::string, std::vector<double>> ppowers;
    std::map<std::string, std::vector<double>> pterms;
    ppowers.insert_or_assign("P1", x);
    std::vector<double> pvalue1, pvalue2, padd_res, pmult_res;

    // 평가
    for(auto step: es.eval_step)
    {
        switch(step.op)
        {
            case 'o':
                if(step.key1[0] != 'C')
                {
                    if(step.key1[0] == 'P')
                        continue;
                    throw std::invalid_argument(std::format("Main::RES: invalid key data {} {}", step.key1, step.key2));
                }
                coeffs.insert_or_assign(step.key1, std::stod(step.key2));
                break;
            case '+':
                if(step.key1[0] == 'T' && step.key2[0] == 'C')
                {
                    padd_res = add_vector_plain(pterms.at(step.key1), coeffs.at(step.key2));
                }
                else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                {
                    padd_res = add_vector(pterms.at(step.key1), pterms.at(step.key2));
                }
                else if(step.key1[0] == 'P' && step.key2[0] == 'C')
                {
                    padd_res = add_vector_plain(ppowers.at(step.key1), coeffs.at(step.key2));
                }
                else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                {
                    padd_res = add_vector(ppowers.at(step.key1), pterms.at(step.key2));
                }
                else
                {
                    throw std::invalid_argument(std::format("Main::ADD: invalid key data {} {}", step.key1, step.key2));
                }
                pterms.insert_or_assign(step.save_key, padd_res);
                break;

            case 'x':
                if(step.key1[0] == 'P' && step.key2[0] == 'P')
                {
                    pmult_res = mult_vector(ppowers.at(step.key1), ppowers.at(step.key2));
                }
                else if(step.key1[0] == 'C' && step.key2[0] == 'P')
                {
                    pmult_res = mult_vector_plain(ppowers.at(step.key2), coeffs.at(step.key1));
                }
                else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                {
                    pmult_res = mult_vector(ppowers.at(step.key1), pterms.at(step.key2));
                }
                else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                {
                    pmult_res = mult_vector(pterms.at(step.key1), pterms.at(step.key2));
                }
                else if(step.key1[0] == 'T' && step.key2[0] == 'P')
                {
                    pmult_res = mult_vector(pterms.at(step.key1), ppowers.at(step.key2));
                }
                else
                {
                    throw std::invalid_argument(std::format("Main::MUL: invalid key data {} {}", step.key1, step.key2));
                }
                switch(step.save_key[0])
                {
                    case 'P':
                        ppowers.insert_or_assign(step.save_key, pmult_res);
                        break;
                    case 'T':
                        pterms.insert_or_assign(step.save_key, pmult_res);
                        break;
                    default:
                        break;
                }
                break;
            default:
                break;
        }
    }
    std::vector<double> pres = pterms.at(es.eval_step.back().save_key);
    return pres;
} */