#pragma once
#include "seal/seal.h"
#include "CKKS_params.h"
#include "evaluate.h"
#include "measure_time.h"
#include "error_bound.h"
#include "polyEval_class.h"
#include "parse.h"
#include "math.h"

#include <format>


inline std::map<int, std::map<int, std::pair<int,int>>> required_iter_HELUT = {
    {35, {
        {2, {6, 1}},
        {3, {8, 1}},
        {4, {10, 2}},
        {5, {11, 2}},
        {6, {13, 3}},
    }},
    {50, {
        {2, {7, 1}},
        {3, {9, 1}},
        {4, {11, 1}},
        {5, {12, 2}},
        {6, {14, 2}},
        {7, {16, 2}},
        {8, {18, 2}},
        {9, {20, 2}},
        {10, {22, 2}},
    }},
};

inline map<int, map<int, int>> required_depth_REMEZ = {
    {35, {
        {2, 99}, {3, 99}, {4, 99}, {5, 99}, {6, 99}
    }},
    {50, {
        {{2, 8}, {3, 9}, {4, 11}, {5, 13}, {6, 14}, {7, 15}, {8, 16}, {9, 19}, {10, 20}}
    }}
};

std::unique_ptr<CKKS_params> set_test_params(int p_num, int s_num, int e_num, int N, size_t hwt, bool test_helut, bool test_remez, int depth=0)
{
    //필요한 깊이 계산
    int helut_level = 0, remez_level = 0;
    pair<int, int> iter_HELUT;
    if(test_helut) 
    {
        iter_HELUT = required_iter_HELUT.at(s_num).at(p_num);
        helut_level = 2 + iter_HELUT.first + (2 * iter_HELUT.second);
    }
    if(test_remez)
    {
        remez_level = depth;
    }
    int req_depth = helut_level > remez_level ? helut_level : remez_level;

    //파라미터 정보 출력
    std::cout << "| Interval p:\t2^" << p_num << std::endl;
    std::cout << "| Precision e:\t2^-" << e_num << std::endl;
    std::cout << "| Scale s:\t2^" << s_num << std::endl;
    std::cout << "| Depth:\t" << req_depth << std::endl;
    std::cout << "------------\n";

    // Modulus chain.
    vector<int> modulus = {60};
    for(int i=0; i<req_depth; i++)
        modulus.push_back(s_num);
    modulus.push_back(60);

    return make_unique<CKKS_params>(modulus, p_num, e_num, s_num, N, hwt, iter_HELUT);
}

std::unique_ptr<CKKS_params> set_test_params(int p_num, int s_num, int e_num, int N, size_t hwt, int req_depth)
{
    //파라미터 정보 출력
    std::cout << "| Interval p:\t2^" << p_num << std::endl;
    std::cout << "| Precision e:\t2^-" << e_num << std::endl;
    std::cout << "| Scale s:\t2^" << s_num << std::endl;
    std::cout << "| Depth:\t" << req_depth << std::endl;
    std::cout << "------------\n";

    // Modulus chain.
    vector<int> modulus = {60};
    for(int i=0; i<req_depth; i++)
        modulus.push_back(s_num);
    modulus.push_back(60);

    return make_unique<CKKS_params>(modulus, p_num, e_num, s_num, N, hwt);
}

std::vector<double> generate_points(std::pair<double, double> interval, int count)
{
    double start = interval.first;
    double end = interval.second;
    std::vector<double> p;

    if (count <= 0) return p;
    if (count == 1) {
        p.push_back(start);
        return p;
    }

    double step = (end - start) / (count - 1);
    for (int i = 0; i < count; ++i) {
        p.push_back(start + i * step);
    }
    return p;
}


// 분해식 정보 기반 가암호문 평가
/* fct::Ciphertext eval_poly_fct(shared_ptr<Decomp> dcmp, fct::Ciphertext& fct, ErrBound& eb)
{
    EvalStep es(dcmp);
    std::map<std::string, fct::Ciphertext> powers;
    std::map<std::string, fct::Ciphertext> coeffs;
    std::map<std::string, fct::Ciphertext> terms;
    powers.insert_or_assign("P1", fct);

    for(auto step: es.eval_step)
    {
        switch(step.op)
        {
            case 'o':
                if(step.key1[0] != 'C')
                {
                    if(step.key1[0] == 'P')
                        continue;
                    throw std::invalid_argument(std::format("Main::FRES: invalid key data {} {}", step.key1, step.key2));
                }
                coeffs.insert_or_assign(step.key1, fct::Ciphertext(std::stod(step.key2), fct.pt.size(), eb.Bc, eb.scale));
                break;
            case '+':
                if(step.key1[0] == 'T' && step.key2[0] == 'C')
                {
                    fct::Ciphertext temp = add(terms.at(step.key1), coeffs.at(step.key2));
                    terms.insert_or_assign(step.save_key, temp);
                }
                else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                {
                    terms.insert_or_assign(step.save_key, add(terms.at(step.key1), terms.at(step.key2)));
                }
                else if(step.key1[0] == 'P' && step.key2[0] == 'C')
                {
                    terms.insert_or_assign(step.save_key, add(powers.at(step.key1), coeffs.at(step.key2)));
                }
                else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                {
                    terms.insert_or_assign(step.save_key, add(powers.at(step.key1), terms.at(step.key2)));
                }
                else
                {
                    throw std::invalid_argument(std::format("Main::FADD: invalid key data {} {}", step.key1, step.key2));
                }
                break;
            case 'x':
                if(step.key1[0] == 'P' && step.key2[0] == 'P')
                {
                    powers.insert_or_assign(step.save_key, mult(powers.at(step.key1), powers.at(step.key2), eb));
                }
                else if(step.key1[0] == 'C' && step.key2[0] == 'P')
                {
                    terms.insert_or_assign(step.save_key, mult(coeffs.at(step.key1), powers.at(step.key2), eb, true));
                }
                else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                {
                    terms.insert_or_assign(step.save_key, mult(powers.at(step.key1), terms.at(step.key2), eb));
                }
                else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                {
                    terms.insert_or_assign(step.save_key, mult(terms.at(step.key1), terms.at(step.key2), eb));
                }
                else if(step.key1[0] == 'T' && step.key2[0] == 'P')
                {
                    terms.insert_or_assign(step.save_key, mult(terms.at(step.key1), powers.at(step.key2), eb));
                }
                else
                {
                    throw std::invalid_argument(std::format("Main::FMUL: invalid key data {} {}", step.key1, step.key2));
                }
                break;
            default:
                break;
        }
    }
    return terms.at(es.eval_step.back().save_key);
} */

// remez - 필요한 Cleanse 횟수 계산
/* void test_cleanse_s(int p_num, int e_num, int s_num)
{
    int pre_num = e_num;
    double scale = pow(2, s_num);
    MR_Info mi = parse_remez_coeff(p_num, e_num, false, false);
    int remez_level = 10;
    for(int i=0; i<mi.coeffs.size(); i++)
        remez_level += evaluate_depth(mi.coeffs[i], true);

    vector<int> modulus = {60};
    for (int i = 0; i < remez_level; i++)
        modulus.push_back(s_num);
    modulus.push_back(60);
    CKKS_params pms(modulus, scale, pow(2, 17), 192);

    vector<double> x;
    for(int i=0; i< int(pow(2, p_num)); i++)
        x.push_back(i);

    Ciphertext ct = pms.encrypt(x);
    vector<vector<double>> coeffs = mi.coeffs;
    vector<double> real_result, ckks_result, coeff;
    real_result = x;
    std::cout << "p, e: " << p_num << ", " << e_num << std::endl;
    std::cout << "precision: " << -e_num << std::endl;

    // 1. coeff
    Ciphertext temp_cleanse;
    for (int en=0; en < coeffs.size(); en++)
    {
        coeff = coeffs[en];
        vector<double> corr_factor = calculate_correction_factor(pms, coeff, ct);
        ct = eval(pms, coeff, corr_factor, ct);

        // 2. Cleanse
        double max_err, pre_err;
        temp_cleanse = ct;
        ckks_result = pms.decode_ctxt(ct);
        max_err = compare_result_log(ckks_result, x.size(), false);
        pre_err = max_err;
        if(max_err < -pre_num)
            std::cout << "| Success(n, s, err): " << en+1 << ", " << 0 << ", " << max_err << std::endl;
        else
            std::cout << "| Fail(n, s, err): " << en+1 << ", " << 0 << ", " << max_err << std::endl;

        for(int i=0; i<5; i++)
        {
            temp_cleanse = cleanse(pms, temp_cleanse);
            ckks_result = pms.decode_ctxt(temp_cleanse);
            max_err = compare_result_log(ckks_result, x.size(), false);
            if ((max_err > pre_err) && (max_err > -pre_num))
                break;
            
            if(max_err <= -pre_num)
                std::cout << "| Success(n, s, err): " << en+1 << ", " << i+1 << ", " << max_err << std::endl;
            else
                std::cout << "| Fail(n, s, err): " << en+1 << ", " << i+1 << ", " << max_err << std::endl;
            pre_err = max_err;
        }
    }
    // std::cout << "##################################" << std::endl;
}
*/

// 시간 측정
/* void test_time(int p_num, int s_num, int e_num, int N_num, bool test_helut, bool test_remez)
{

    if(test_remez)
        test_remez_time(p_num, e_num, s_num, N_num, true);
    if(test_helut)
        test_helut_time(p_num, e_num, s_num, N_num, true);
    std::cout << "##########################" << std::endl;
} */

// err bound 측정
/* void test_err(EB eb, vector<double> coeff, int p_num, int e_num, int s_num, int N_num, int hwt)
{
    // 기본 파라미터 설정
    size_t poly_modulus_degree = pow(2, N_num);

    double p = pow(2, p_num);
    double e = pow(2, -e_num);
    double s = pow(2, s_num);

    int remez_level = evaluate_depth(coeff, true);

    // Modulus chain.
    vector<int> modulus = {60};
    for(int i=0; i<remez_level; i++)
        modulus.push_back(s_num);
    modulus.push_back(60);

    CKKS_params pms(modulus, s, poly_modulus_degree, hwt);

    vector<double> x;
    for(int i=0; i<pow(2, p_num); i++)
        x.push_back(i); 
    Ciphertext ct = pms.encrypt(x);
    vector<double> corr_factor = calculate_correction_factor(pms, coeff, ct);
    ct = eval(pms, coeff, corr_factor, ct);

    vector<double> real_result = evaluate_function(coeff, x);
    vector<double> ckks_result = pms.decode_ctxt(ct);
    double err_bound, err;

    //결과 비교
    err_bound = eb.cal_bound(pow(2, p_num)-1, coeff)/pms.scale;
    cout << "Bc: " << eb.Bc << endl;
    cout << "Bs: " << eb.Bs << endl;
    for(int i=0; i<real_result.size(); i++)
    {
        err = abs(real_result[i] - ckks_result[i]);
        if(err > err_bound)
            cout << "ERROR OUTSIDE BOUND!!!!" << endl;
        
        cout << real_result[i] << "\t\t\t\t" << ckks_result[i] << endl;
        cout << err_bound << "\t\t\t" << err << "\n" << endl;
    }
    cout << "------------\n";
}
*/