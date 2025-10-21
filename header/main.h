#pragma once
#include "seal/seal.h"
#include "CKKS_params.h"
#include "parse.h"
#include "evaluate.h"
#include "measure_time.h"
#include "tbench.h"
#include "error_bound.h"
#include <cmath>

using namespace std;
using namespace seal;


void test_base(int p_num, int s_num, int e_num, int N, size_t hwt, bool test_helut, bool test_remez, bool mid_print, bool res_print)
{
    // 기본 파라미터 설정
    size_t poly_modulus_degree = pow(2, N);
    std::chrono::_V2::system_clock::time_point start_time, end_time;
    int req_depth;

    double p = pow(2, p_num);
    double e = pow(2, -e_num);
    double s = pow(2, s_num);

    HELUT_Info hi;
    MR_Info mi;

    //필요한 깊이 계산
    int helut_level = 0, remez_level = 0;
    if(test_helut)
    {
        hi = parse_helut(p_num, e_num);
        helut_level = 2 + hi.r + 2 * hi.s;
    }
    if(test_remez)
    {
        mi = parse_remez_coeff(p_num, e_num, true, false);
        for(int i=0; i<mi.n - mi.s; i++)
        {
            vector<double> coeff = mi.coeffs[i];
            remez_level += evaluate_depth(coeff, true);
        }           
        remez_level += (2 * mi.s);
    }
    req_depth = helut_level > remez_level ? helut_level : remez_level;

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

    CKKS_params pms(modulus, s, poly_modulus_degree, hwt);

    // HELUT 평가
    if(test_helut)
    {
        std::cout << "Evaluating HELUT" << std::endl;
        start_time = cur_time();
        evaluate_helut(pms, hi, e, p, mid_print, res_print);
        end_time = cur_time();
        calculate_time(start_time, end_time);
        std::cout << "#########################################################" << std::endl;
    }

    // Remez 평가
    if(test_remez)
    {
        std::cout << "Evaluating Remez" << std::endl;
        start_time = cur_time();
        evaluate_multi_remez(pms, mi, e, p, mid_print, res_print);
        end_time = cur_time();
        calculate_time(start_time, end_time);
    }
}



// remez - 필요한 Cleanse 횟수 계산
/*
void test_cleanse_s(int p_num, int e_num, int s_num)
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
void test_time(int p_num, int s_num, int e_num, int N_num, bool test_helut, bool test_remez)
{

    if(test_remez)
        test_remez_time(p_num, e_num, s_num, N_num, true);
    if(test_helut)
        test_helut_time(p_num, e_num, s_num, N_num, true);
    std::cout << "##########################" << std::endl;
}


// err bound 측정
/*
void test_err(EB eb, vector<double> coeff, int p_num, int e_num, int s_num, int N_num, int hwt)
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