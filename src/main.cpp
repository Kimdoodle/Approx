#include "../include/main.h"

const int p_num = 10;
const int e_num = 30;
const int s_num = 50;
const int N_num = 17;
const size_t hwt = 192;
const double sigma = 3.1;
const double scale = pow(2, s_num);
const int N = pow(2, N_num);
const bool step_debug = false;
const bool coeff_debug = false;
const bool res_debug = true;

const std::vector<double> cl_coeff = {0, 0, 3.0, -2.0};

// 특정 다항식의 오차범위가 정확히 계산되는지 확인
/* int main()
{
    // 기본 파라미터, 분해식정보 불러오기
    ErrBound eb(sigma, N, hwt, s_num);

    std::vector<double> coeff = {0, 0, 3, -2};

    vector<std::shared_ptr<Decomp>> dcmps;
    json decomp_cache = load_decomp_cache("data/decomp_cache.json");
    std::unique_ptr<CKKS_params> pms = set_test_params(p_num, s_num, e_num, N_num, hwt, ceil(log2(int(coeff.size()))));
    
    // 기본 평문 데이터 설정
    int n_value = pow(2, N_num);
    int p_value = pow(2, p_num);
    std::vector<double> x;
    x.reserve(n_value/2);
    int copy_count = x.capacity() / p_value;
    for(int i=0; i<p_value; i++)
        for(int j=0; j<copy_count; j++)
            x.push_back(static_cast<double>(i));
    seal::Ciphertext ctx = pms->encrypt(x);

    const int sample_count = copy_count;
    std::vector<double> x2;
    x2.reserve(sample_count * p_value);
    std::vector<double> points;
    const double bc = eb.Bc/eb.scale;
    for(int i=0; i<p_value; i++)
    {
        points = generate_points({i-bc, i+bc}, sample_count);
        x2.insert(x2.end(), points.begin(), points.end());
    }
    fct::Ciphertext fctx(x2, eb.Bc, eb.scale);

    // 다항식 평가
    std::tuple<std::vector<double>, seal::Ciphertext, fct::Ciphertext> res = {x, ctx, fctx};
    string key = get_poly_type_key(coeff);
    shared_ptr<Decomp> dcmp = reconstruct_decomp_from_cache(decomp_cache.at(key), coeff);
    dcmps.push_back(dcmp);
    // EvalStep es(dcmp);

    // pt, ct, fct 한번에 평가
    res = eval_poly_ct(*pms, dcmp, res, eb, copy_count, true);

    return 0;
} */

int main()
{
    // 기본 파라미터, 분해식정보 불러오기
    ErrBound eb(sigma, N, hwt, s_num);
    vector<vector<double>> coeffs = parse_remez_coeff(p_num, e_num);
    int depth = calcuate_remez_depth(coeffs);
    vector<std::shared_ptr<Decomp>> dcmps;
    json decomp_cache = load_decomp_cache("data/decomp_cache.json");
    std::unique_ptr<CKKS_params> pms = set_test_params(p_num, s_num, e_num, N_num, hwt, false, true, depth);
    
    // 기본 평문 데이터 설정
    // 각 x값을 최대로 삽입.
    int n_value = pow(2, N_num);
    int p_value = pow(2, p_num);
    std::vector<double> x;
    x.reserve(n_value / 2);
    int copy_count = x.capacity();
    for (int i = 0; i < p_value; ++i)
    {
        std::fill(x.begin(), x.end(), static_cast<double>(i));
    
        seal::Ciphertext ctx = pms->encrypt(x);
        fct::Ciphertext fctx(x, eb.Bc, eb.scale);

        // 각 다항식에 대해 평가
        std::vector<double> res_vector = x;
        seal::Ciphertext res_ct = ctx;
        fct::Ciphertext res_fct = fctx;
        std::tuple<std::vector<double>, seal::Ciphertext, fct::Ciphertext> res = {x, ctx, fctx};
        for(vector<double> coeff: coeffs)
        {
            if(coeff == cl_coeff)
            {
                if(step_debug || coeff_debug)
                    std::cout << std::format("Polynomial\t{}\n", vec_to_str(coeff));
                res = eval_poly_cl(*pms, res, eb, copy_count, step_debug, coeff_debug);
            }
            else
            {
                string key = get_poly_type_key(coeff);
                shared_ptr<Decomp> dcmp = reconstruct_decomp_from_cache(decomp_cache.at(key), coeff);
                dcmps.push_back(dcmp);
                if(step_debug || coeff_debug)
                    std::cout << std::format("Polynomial\t{}\n", dcmp->restore_dcmp());

                // pt, ct, fct 한번에 평가
                res = eval_poly_ct(*pms, dcmp, res, eb, copy_count, step_debug, coeff_debug);
            }
            if(step_debug || coeff_debug)
                std::cout << "#####################################" << std::endl;
        }
        if(res_debug)
        {
            std::cout << std::format("Results of x={}\n", i);
            compare_precision(res, *pms, copy_count);
        }
    }
}


/* int main()
{
    bool test_helut = true;
    bool test_remez = false;
    bool test_time = false;
    ErrBound eb(sigma, N, hwt, s_num);

    std::chrono::_V2::system_clock::time_point start_time, end_time;

    // 다항식 데이터 불러오기
    if(test_remez)
    {
        vector<vector<double>> coeffs = parse_remez_coeff(p_num, e_num);
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
    }

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
        // CKKS 파라미터 설정
        std::unique_ptr<CKKS_params> pms = set_test_params(p_num, s_num, e_num, N_num, hwt, test_helut, false);
        std::cout << "Evaluating HELUT" << std::endl;
        start_time = cur_time();
        evaluate_helut(*pms, x, test_time);
        end_time = cur_time();
        if(test_time)
            calculate_time(start_time, end_time);
        std::cout << "#########################################################" << std::endl;
    }

    // Remez 평가
    if(test_remez)
    {
        // CKKS 파라미터 설정
        std::unique_ptr<CKKS_params> pms = set_test_params(p_num, s_num, e_num, N_num, hwt, false, true);
        std::cout << "Evaluating REMEZ" << std::endl;
        start_time = cur_time();
        // evaluate_multi_remez(*pms, x, dcmps, test_time);
        end_time = cur_time();
        if(test_time)
            calculate_time(start_time, end_time);
    }
}
 */