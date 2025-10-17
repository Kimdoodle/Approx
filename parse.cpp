// parse.cpp
#include "header/parse.h"
#include <fstream>
#include <regex>
#include <stdexcept>
#include <iostream>
#include <map>
#include <limits>

// Approx_Data::Approx_Data(int p_num, int s_num)
// {
//     // 1. HELUT
//     string filename_helut = "doc/helut.txt";
//     this->hi = parse_helut(p_num, s_num);
//     // 3. Remez
//     string filename_remez = "doc/coeff_" + to_string(p_num) + "_" + to_string(nextPowerOfTwo(hi.e)) + ".txt";
//     this->mi = parse_remez_coeff(filename_remez);
// }
static double stod_safe(const std::string& s) {
    // 필요시 안전 변환용 헬퍼
    return std::stod(s);
}

HELUT_Info parse_helut(int p_num, int e_num) {
    string filename = "doc/helut.txt";
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("파일을 열 수 없습니다: " + filename);

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string token;
        int p = 0, d = 0, e = 0, r = 0, u = 0, s = 0;

        while (iss >> token) {
            if (token.rfind("p=", 0) == 0)      p = std::stoi(token.substr(2));
            else if (token.rfind("d=", 0) == 0) d = std::stoi(token.substr(2));
            else if (token.rfind("e=", 0) == 0) e = std::stoi(token.substr(2));
            else if (token.rfind("r=", 0) == 0) r = std::stoi(token.substr(2));
            else if (token.rfind("u=", 0) == 0) u = std::stoi(token.substr(2));
            else if (token.rfind("s=", 0) == 0) s = std::stoi(token.substr(2));
        }

        if (p == p_num && e == e_num) {
            return HELUT_Info{r, s};
        }
    }

    throw std::runtime_error("일치하는 (p, e) 데이터를 찾을 수 없습니다.");
}

MR_Info parse_remez_coeff(int p_num, int e_num, bool parse_s, bool print)
{
    int s_val = 0;  // "-2x^3+3x^2"가 포함된 줄의 개수
    int n_val = 0;  // 전체 줄 개수

    // 1) coeff 파일명을 (p,e)로 구성
    const std::string filename = "doc/coeff_" + std::to_string(p_num) + "_" + std::to_string(e_num) + ".txt";

    // 2) 계수 파일 파싱
    std::ifstream infile(filename);
    if (!infile)
        throw std::runtime_error("파일을 열 수 없습니다: " + filename);

    MR_Info mi;
    std::vector<std::vector<double>> result;
    std::string line;

    // 계수 추출 패턴: (계수)x^(지수)
    const std::regex term_re(
        R"(([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?))" // coeff (그룹1)
        R"(x\^)"                                              // x^
        R"((\d+))"                                            // exponent (그룹2)
    );

    // s 계산용 패턴: "-2x^3+3x^2" (공백 가변 허용)
    const std::regex s_re(R"(-\s*2\s*x\^3\s*\+\s*3\s*x\^2)");

    // 파일 순회
    while (std::getline(infile, line)) {
        ++n_val;  // 줄 개수 증가

        // parse_s가 true인 경우 "-2x^3+3x^2" 포함 여부 확인
        if (parse_s && std::regex_search(line, s_re)) {
            ++s_val;
        }

        // (exp -> coeff) 맵 구성
        std::map<int, double> exp_to_coeff;
        int max_exp = -1;

        for (std::sregex_iterator it(line.begin(), line.end(), term_re);
             it != std::sregex_iterator(); ++it)
        {
            const std::string coeff_str = (*it)[1].str();
            const std::string exp_str   = (*it)[2].str();

            double c = stod_safe(coeff_str);
            int    e = std::stoi(exp_str);

            exp_to_coeff[e] = c;
            if (e > max_exp) max_exp = e;
        }

        // 항이 하나도 없으면 빈 벡터
        if (max_exp < 0) {
            result.emplace_back();
            continue;
        }

        // 누락 차수 0.0으로 채움
        std::vector<double> coeffs(static_cast<size_t>(max_exp + 1), 0.0);
        for (const auto& kv : exp_to_coeff) {
            coeffs[static_cast<size_t>(kv.first)] = kv.second;
        }

        result.push_back(std::move(coeffs));
    }

    // 3) 출력 옵션
    if (print) {
        std::cout << "[parse_remez_coeff] p=" << p_num
                  << ", e=" << e_num
                  << ", n=" << n_val
                  << ", s=" << s_val << "\n";
        for (size_t i = 0; i < result.size(); ++i) {
            std::cout << "Line " << i << " coeffs (k: coeff for x^k):\n  ";
            const auto& v = result[i];
            for (size_t k = 0; k < v.size(); ++k) {
                std::cout << "x^" << k << ":" << v[k];
                if (k + 1 != v.size()) std::cout << ", ";
            }
            std::cout << "\n";
        }
    }

    // 4) MR_Info 채우기
    mi.coeffs = std::move(result);
    mi.n = n_val;   // 읽은 줄 수
    mi.s = s_val;   // "-2x^3+3x^2"가 포함된 줄 수

    return mi;
}

// MR_Info parse_remez_coeff(int p_num, int e_num, bool parse_s, bool print)
// {
//     int s_val = 0;
//     int n_val = 0;
//     if(parse_s)
//     {
//         // 0) remez.txt에서 (p,e)에 대응하는 s 값 읽기
//         const std::string meta_path = "doc/remez.txt";
//         std::ifstream meta_in(meta_path);
//         if (!meta_in)
//             throw std::runtime_error("파일을 열 수 없습니다: " + meta_path);

//         // 형식: "p, e, s: _, _, _"  (밑줄은 정수)
//         // 공백은 가변적이라고 가정하고 정규식으로 파싱
//         const std::regex re_meta(
//             R"(p\s*,\s*e\s*,\s*n\s*,\s*s\s*:\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(-?\d+))"
//         );
//         bool found_meta = false;
//         std::string mline;
//         while (std::getline(meta_in, mline)) {
//             std::smatch m;
//             if (std::regex_search(mline, m, re_meta)) {
//                 const int p = std::stoi(m[1].str());
//                 const int e = std::stoi(m[2].str());
//                 const int n = std::stoi(m[3].str());
//                 const int s = std::stoi(m[4].str());
//                 if (p == p_num && e == e_num) {
//                     n_val = n;
//                     s_val = s;
//                     found_meta = true;
//                     break;
//                 }
//             }
//         }
//         if (!found_meta) {
//             throw std::runtime_error(
//                 "remez.txt에서 주어진 (p,e)=(" + std::to_string(p_num) + "," + std::to_string(e_num) + ") 에 해당하는 s를 찾을 수 없습니다."
//             );
//         }
//     }

//     // 1) coeff 파일명을 (p,e)로 구성
//     const std::string filename = "doc/coeff_" + std::to_string(p_num) + "_" + std::to_string(e_num) + ".txt";

//     // 2) 계수 파일 파싱 (기존 로직)
//     std::ifstream infile(filename);
//     if (!infile)
//         throw std::runtime_error("파일을 열 수 없습니다: " + filename);

//     MR_Info mi;
//     std::vector<std::vector<double>> result;
//     std::string line;

//     // 패턴: (계수)x^(지수)
//     // 예: -1.2315e-05x^0, +3.9924x^2, 7x^10
//     const std::regex term_re(
//         R"(([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?))" // coeff (그룹1)
//         R"(x\^)"                                              // x^
//         R"((\d+))"                                            // exponent (그룹2)
//     );

//     while (std::getline(infile, line)) {
//         // (exp -> coeff) 맵 구성
//         std::map<int, double> exp_to_coeff;
//         int max_exp = -1;

//         for (std::sregex_iterator it(line.begin(), line.end(), term_re);
//              it != std::sregex_iterator(); ++it)
//         {
//             const std::string coeff_str = (*it)[1].str();
//             const std::string exp_str   = (*it)[2].str();

//             double c = stod_safe(coeff_str);
//             int    e = std::stoi(exp_str);

//             // 동일 차수 중복 시 마지막 항으로 갱신
//             exp_to_coeff[e] = c;
//             if (e > max_exp) max_exp = e;
//         }

//         // 항이 하나도 없으면 빈 벡터(빈 줄 처리)
//         if (max_exp < 0) {
//             result.emplace_back();
//             continue;
//         }

//         // 누락 차수 0.0으로 채운 계수 벡터 구성
//         std::vector<double> coeffs(static_cast<size_t>(max_exp + 1), 0.0);
//         for (const auto& kv : exp_to_coeff) {
//             coeffs[static_cast<size_t>(kv.first)] = kv.second;
//         }

//         result.push_back(std::move(coeffs));
//     }

//     // 3) 출력 옵션
//     if (print) {
//         std::cout << "[parse_remez_coeff] p=" << p_num
//                   << ", e=" << e_num
//                   << ", s=" << s_val << "\n";
//         for (size_t i = 0; i < result.size(); ++i) {
//             std::cout << "Line " << i << " coeffs (k: coeff for x^k):\n  ";
//             const auto& v = result[i];
//             for (size_t k = 0; k < v.size(); ++k) {
//                 std::cout << "x^" << k << ":" << v[k];
//                 if (k + 1 != v.size()) std::cout << ", ";
//             }
//             std::cout << "\n";
//         }
//     }

//     // 4) MR_Info 채우기
//     mi.coeffs = std::move(result);
//     mi.n = n_val;
//     mi.s = s_val;

//     return mi;
// }

// MR_Info parse_remez_coeff(int p_num, int e_num, bool print) { 
//     string filename_remez = "doc/coeff_" + to_string(p_num) + "_" + to_string(e_num) + ".txt";
//     string filename2 = "doc/remez.txt";
//     std::ifstream infile(filename_remez); 
//     if (!infile) 
//         throw std::runtime_error("파일을 열 수 없습니다: " + filename_remez); 
//     MR_Info mi; std::vector<std::vector<double>> result; 
//     std::string line; 
    
//     // (계수)x^(지수) 패턴 (부호/소수/지수표기 허용) 
//     // 예: -1.2315e-05x^0, +3.9924x^2, 7x^10 
//     const std::regex term_re( 
//         R"(([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?))" // coeff (그룹1) 
//         R"(x\^)" // x^ 
//         R"((\d+))" // exponent (그룹2) 
//     ); 
//     while (std::getline(infile, line)) 
//     { 
//         // 1) 한 줄에서 (exp -> coeff)로 모으기 
//         std::map<int, double> exp_to_coeff; 
//         int max_exp = -1; 
//         for (std::sregex_iterator it(line.begin(), line.end(), term_re); it != std::sregex_iterator(); ++it) 
//         { 
//             const std::string coeff_str = (*it)[1].str(); 
//             const std::string exp_str = (*it)[2].str(); 
//             double c = stod_safe(coeff_str); 
//             int e = std::stoi(exp_str); 
            
//             // 동일 차수 중복 시 마지막 항으로 갱신 
//             exp_to_coeff[e] = c; 
//             if (e > max_exp) 
//                 max_exp = e; 
//         } 
        
//         // 2) 항이 하나도 없으면 빈 벡터를 추가(또는 스킵하려면 continue) 
//         if (max_exp < 0) 
//         { 
//             result.emplace_back(); 
//             // 빈 줄 처리 
//             continue; 
//         } 
        
//         // 3) 누락 차수 0.0으로 채운 계수 벡터 구성 
//         std::vector<double> coeffs(static_cast<size_t>(max_exp + 1), 0.0); 
//         for (const auto& kv : exp_to_coeff) 
//         { 
//             coeffs[static_cast<size_t>(kv.first)] = kv.second; 
//         } 
//         result.push_back(std::move(coeffs)); 
//     } 
    
//     if (print) 
//     { 
//         for (size_t i = 0; i < result.size(); ++i) 
//         { 
//             std::cout << "Line " << i << " coeffs (k: coeff for x^k):\n "; 
//             const auto& v = result[i]; 
//             for (size_t k = 0; k < v.size(); ++k) 
//             { 
//                 std::cout << "x^" << k << ":" << v[k]; 
//                 if (k + 1 != v.size()) std::cout << ", "; 
//             } 
//             std::cout << "\n"; 
//         } 
//     } 
//     mi.coeffs = result; 
//     return mi; 
// }


unsigned int nextPowerOfTwo(unsigned int n) {
    if (n == 0) return 1;

    // n이 이미 2의 거듭제곱이면 그대로 반환
    if ((n & (n - 1)) == 0) return n;

    // 비트를 오른쪽으로 계속 OR 하여 최상위 비트 아래를 모두 1로 만든 후 +1
    unsigned int p = 1;
    while (p < n) {
        p <<= 1;
    }
    return p;
}