#include "../include/main.h"
#include <iostream>
#include <fstream>

// nohup ./build/APPROX > /dev/null 2>&1 &

const int e_num = 30;
const int s_num = 50;

const double scale = pow(2, s_num);

const bool step_debug = false;
const bool coeff_debug = false;
const bool res_debug = true;

// 다항식의 평가순서 계산
/* int main()
{
    std::vector<double> coeff = {1.2, 3.4, 5.6, 7.8, 9.0, 1.01, 8.45, 2.8, 4.9};
    std::vector<std::shared_ptr<Decomp>> dcmps;
    json decomp_cache = load_decomp_cache("data/decomp_cache.json");
    string key = get_poly_type_key(coeff);
    std::cout << key << std::endl;
    shared_ptr<Decomp> dcmp = reconstruct_decomp_from_cache(decomp_cache.at(key), coeff);
    dcmps.push_back(dcmp);
    EvalStep es(dcmp);
    es.print_step();
} */

// 특정 다항식의 오차범위가 정확히 계산되는지 확인
/* int main()
{
    // 기본 파라미터, 분해식정보 불러오기
    int p_num = 8;
    int sample_num = std::pow(2, 14);

    ErrBound eb(sigma, N, hwt, s_num);

    std::vector<double> coeff = parse_remez_coeff(p_num, e_num, "depth")[0];
    std::vector<std::shared_ptr<Decomp>> dcmps;
    json decomp_cache = load_decomp_cache("data/decomp_cache.json");
    EvalStep es = EvalStep(reconstruct_decomp_from_cache(decomp_cache.at(get_poly_type_key(coeff)), coeff));
    
    double start = 0.0;
    double p_value = pow(2, p_num);
    double err = eb.Bc / eb.scale;
    std::vector<double> max_interval_size(p_value);
    for(int i=0; i<=p_value; i++)
    {
        //fct::Ciphertext
        fct::Ciphertext ct = fct::Ciphertext(i-err, i+err, sample_num, 0.0, eb.scale, true);
        ct = evaluate(ct, es, eb);

        auto [res_fct_max, res_fct_min] = fct::dec(ct);
        for(int i=0; i<res_fct_max.size(); i++)
        {
            assert(res_fct_max[i] > res_fct_min[i]), "Invalid interval: max value is smaller than min value!";
        }
        // debug - 구간 크기 측정
        std::vector<double> slot_size = pt::sub(res_fct_max, res_fct_min);
        double fct_min = *std::min_element(res_fct_min.begin(), res_fct_min.end());;
        double fct_max = *std::max_element(res_fct_max.begin(), res_fct_max.end());
        max_interval_size[i] = *std::max_element(slot_size.begin(), slot_size.end());

        // seal::Ciphertext
        std::unique_ptr<CKKS_params> pms = set_test_params(p_num, s_num, e_num, N_num, hwt, ceil(log2(int(coeff.size()))));
        int n_value = pow(2, N_num);
        std::vector<double> ptx(n_value/2, (double)i);
        seal::Ciphertext ctx = pms->encrypt(ptx);
        ctx = ct::evaluate(ctx, es, *pms);
        std::vector<double> res_pt = pms->decode_ctxt(ctx);
        double ct_min = *std::min_element(res_pt.begin(), res_pt.end());
        double ct_max = *std::max_element(res_pt.begin(), res_pt.end());

        if((fct_min >= ct_min))
        {
            std::cout << std::format("Min boundary error in x={}!\n{} > {}, margin={}\n", i, clearNum(fct_min), clearNum(ct_min), clearNum(abs(fct_min - ct_min)));
        }
        if((fct_max <= ct_max))
        {
            std::cout << std::format("Max boundary error in x={}!\n{} < {}, margin={}\n", i, clearNum(fct_max), clearNum(ct_max), clearNum(abs(fct_max - ct_max)));
        }
    }
    auto max_it = std::max_element(max_interval_size.begin(), max_interval_size.end());
    int max_idx = static_cast<int>(std::distance(max_interval_size.begin(), max_it));
    auto min_it = std::min_element(max_interval_size.begin(), max_interval_size.end());
    int min_idx = static_cast<int>(std::distance(max_interval_size.begin(), min_it));

    std::cout << std::format("max interval size: {} in x={}", clearNum(*max_it), max_idx) << std::endl;
    std::cout << std::format("min interval size: {} in x={}", clearNum(*min_it), min_idx) << std::endl;
    return 0;
} */

int main()
{
    // Margin test
    margin_test(9, 9, true, s_num, e_num, true);
    return 0;
}

// Running time test
/* int main()
{
    std::ofstream ofs(std::format("runtime_test_{}.txt", e_num));
    std::streambuf* old_buf = std::cout.rdbuf(ofs.rdbuf());
    std::streambuf* old_cerr_buf = std::cerr.rdbuf(ofs.rdbuf());

    int test_round = 100;
    std::string critieria = "depth";
    auto total_time_us = std::chrono::microseconds::zero();
    bool test_helut = false;
    bool test_remez = true;
    bool test_time = true;
    std::chrono::_V2::system_clock::time_point start_time, end_time;

    for(int p_num=10; p_num<=10; p_num++)
    {
        std::cout << std::format("TEST P={}", p_num) << std::endl;
        // 평문 데이터 구성 - 동일한 데이터를 가득 채우기
        int n_value = pow(2, N_num);    
        int p_value = pow(2, p_num);
        std::vector<double> x;
        x.reserve(n_value/2);
        int copy_count = x.capacity() / p_value;
        for(int i=0; i<p_value; i++)
            for(int j=0; j<copy_count; j++)
                x.push_back(static_cast<double>(i));

        // HELUT 평가
        if(test_helut)
        {
            total_time_us = std::chrono::microseconds::zero();
            for(int i=0; i<test_round; i++)
            {
                std::cout << std::format("HELUT evaluation - round {} / {}\t", i+1, test_round);
                std::unique_ptr<CKKS_params> pms = set_test_params(p_num, s_num, e_num, N_num, hwt, test_helut, false, false);
                start_time = cur_time();
                evaluate_helut(*pms, x, test_time);
                end_time = cur_time();
                total_time_us += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                calculate_time(start_time, end_time);
            }
            std::cout << std::format("Average HELUT evaluation time: {} s\n", total_time_us.count() / test_round / 1000000.0);
        }

        // Remez 평가
        if(test_remez)
        {
            vector<vector<double>> coeffs = parse_remez_coeff(p_num, e_num, critieria);
            vector<std::shared_ptr<Decomp>> dcmps;

            //분해식 데이터 파싱
            string dcmp_cache_filename = "data/decomp_cache.json";
            json decomp_cache = load_decomp_cache(dcmp_cache_filename);
            if(decomp_cache != nullptr)
                cout << "Cache loaded successfully." << endl << endl;
            else
                return 0;
            for(vector<double> coeff: coeffs)
            {
                string key = get_poly_type_key(coeff);
                shared_ptr<Decomp> result = reconstruct_decomp_from_cache(decomp_cache.at(key), coeff);
                dcmps.push_back(result);
            }

            total_time_us = std::chrono::microseconds::zero();
            for(int i=0; i<test_round; i++)
            {
                int depth = calcuate_remez_depth(coeffs);
                std::unique_ptr<CKKS_params> pms = set_test_params(p_num, s_num, e_num, N_num, hwt, depth);
                // std::unique_ptr<CKKS_params> pms = set_test_params(p_num, s_num, e_num, N_num, hwt, false, true, false);
                seal::Ciphertext temp = pms->encrypt(x);
                std::cout << std::format("REMEZ evaluation - round {} / {}\t", i+1, test_round);
                start_time = cur_time();
                // evaluate_multi_remez(*pms, x, dcmps, test_time);
                for(auto dcmp: dcmps)
                {
                    temp = eval_poly_ct(*pms, dcmp, temp);
                }
                end_time = cur_time();
                total_time_us += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                calculate_time(start_time, end_time);
            }
            std::cout << std::format("Average REMEZ evaluation time: {} s\n", total_time_us.count() / test_round / 1000000.0);
        }
        std::cout << std::endl;
    }
} */
