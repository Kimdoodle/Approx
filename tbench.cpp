// 테스트용 코드
#include "header/tbench.h"
using namespace std;

// 1. HELUT, REMEZ 실행시간 비교

// HELUT 시간 측정
void test_helut_time(int p_num, int e_num, int s_num, bool print)
{
    int iter = 10;
    chrono::_V2::system_clock::time_point start_time, end_time;
    chrono::milliseconds duration;
    double total_seconds = 0.0;
    double avg_seconds, seconds;
    double s = pow(2, s_num);
    int N = 17;
    size_t poly_modulus_degree = pow(2, N);
    size_t hwt = 192;
    HELUT_Info hi = parse_helut(p_num, e_num);
    
    // modulus chain
    vector<int> modulus = {60};
    for(int i=0; i<2 + hi.r + 2*hi.s; i++)
        modulus.push_back(s_num);
    modulus.push_back(60);

    for(int i=0; i<iter; i++)
    {
        poly_modulus_degree = pow(2, N);
        CKKS_params pms(modulus, s, poly_modulus_degree, hwt);
        start_time = cur_time();
        evaluate_helut(pms, hi, pow(2, -e_num), pow(2, p_num), false, false);
        end_time = cur_time();
        duration = duration_cast<milliseconds>(end_time - start_time);
        seconds = duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
        total_seconds += seconds;
    }
    avg_seconds = total_seconds / iter;
    if(print)
    {
        cout << "| Interval p:\t2^" << p_num << endl;
        cout << "| Precision e:\t2^-" << e_num << endl;
        cout << "| Scale s:\t2^" << 50 << endl;
        cout << "| Depth:\t" << 2 + hi.r + 2*hi.s << endl;
        cout << "| Modulus:\t2^" << N << endl;
        cout << "| Average time:\t" << avg_seconds << "s" << endl;
        cout << "------------\n";
    }
}

// Remez 시간 증가량 측정
void test_remez_time(int p_num, int e_num, int s_num, bool print)
{
    int iter = 10;
    chrono::_V2::system_clock::time_point start_time, end_time;
    chrono::milliseconds duration;
    double total_seconds = 0.0;
    double avg_seconds, seconds;
    double s = pow(2, s_num);
    int pmd;
    size_t poly_modulus_degree;

    MR_Info mi = parse_remez_coeff(p_num, e_num, true);

    // depth count
    int remez_level = 0;
    for(int i=0; i<mi.n; i++)
        remez_level += evaluate_depth(mi.coeffs[i], true);
    // for(vector<double> coeff: mi.coeffs)
    //     remez_level += evaluate_depth(coeff, true);
    remez_level += 2 * mi.s;

    // modulus chain
    vector<int> modulus = {60};
    for(int i=0; i<remez_level; i++)
        modulus.push_back(50);
    modulus.push_back(60);

    pmd = 17;
    poly_modulus_degree = pow(2, pmd);
    size_t hwt = 192;
    CKKS_params pms(modulus, s, poly_modulus_degree, hwt);
    for(int i=0; i<iter; i++)
    {
        start_time = cur_time();
        evaluate_multi_remez(pms, mi, pow(2, -e_num), pow(2, p_num), false, false);
        end_time = cur_time();
        duration = duration_cast<milliseconds>(end_time - start_time);
        seconds = duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
        total_seconds += seconds;
        // cout << seconds << endl;
    }
    avg_seconds = total_seconds / iter;
    if(print)
    {
        cout << "| Interval p:\t2^" << p_num << endl;
        cout << "| Precision e:\t2^-" << e_num << endl;
        cout << "| Scale s:\t2^" << 50 << endl;
        cout << "| Depth:\t" << remez_level << endl;
        cout << "| Modulus:\t2^" << pmd << endl;
        cout << "| Average time:\t" << avg_seconds << "s" << endl;
        cout << "------------\n";
    }
}


int find_best_pmd(int s, int depth)
{
    int max_bit_count;
    int i = 10;
    while(i<=15)
    {
        size_t pmd = pow(2, i);
        max_bit_count = CoeffModulus::MaxBitCount(pmd, sec_level_type::tc128);
        if((120 + s*depth) > max_bit_count)
            ++i;
        else
            return i;
    }
    // i=16, 17
    if((120 + s*depth) < 881*2)
        return 16;
    else
        return 17;
}