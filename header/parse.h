#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <regex>

using namespace std;

struct HELUT_Info
{
    int r, s;
};

struct SGN_Info
{
    vector<double> coeff; // 사용할 sgn함수 설정
    int iter; // 반복 횟수
};

struct MR_Info 
{
    vector<vector<double>> coeffs;
    int n;
    int s;
};


class Approx_Data
{
public:
    HELUT_Info hi;
    SGN_Info si;
    MR_Info mi;
    
    // p, s 정보 입력 시 대응하는 HELUT정보 추출
    Approx_Data(int p_num, int s_num);
};

HELUT_Info parse_helut(int p_num, int e_num);
MR_Info parse_remez_coeff(int p_num, int e_num, bool parse_s, bool print = false);
unsigned int nextPowerOfTwo(unsigned int n);
