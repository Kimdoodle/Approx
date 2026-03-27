// decomp.h
#pragma once

#include <vector>
#include <string>
#include <memory>
#include <set>
#include <deque>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "json.hpp"
#include <format>
#define ZERO_CRITERIA 1e-16

using json = nlohmann::json;

// 계산복잡도
class Complexity {
public:
    int depth = 0;
    int cmult = 0;
    int pmult = 0;
    int add = 0;
    
    Complexity();

    void insert_value(int depth, int cmult, int pmult, int add);
};

// 다항식 정보
class Poly {
public:
    std::vector<double> coeff;
    int deg;
    Complexity complexity;

    Poly(std::vector<double> coeff);

    // seperate - f(x) -> (x^i)p(x)+q(x)로 분리.
    std::pair<Poly, Poly> seperate_poly(int i, bool multA);
};

// x^i 정보
class XI {
public:
    bool multA;
    int n;
    int add_count;
    double coeff;
    std::vector<std::pair<int, int>> route;
    std::set<int> made_powers;
    int depth, pmult;

    XI(bool multA=false, int n=0);

    void add_routes(std::vector<std::pair<int, int>> route, std::set<int> made_powers);
};

// 분해식 클래스
class Decomp {
public:
    std::vector<double> coeff;
    // Complexity comp;
    XI xi;
    std::shared_ptr<Decomp> dcmp_p, dcmp_q;
    std::set<int> made_powers;
    std::vector<std::pair<int, int>> total_routes = {};

    Decomp(std::vector<double>& c, Complexity comp, XI xi=XI());

    void update(XI xi, std::shared_ptr<Decomp> dcmp_p, std::shared_ptr<Decomp> dcmp_q);
    std::string restore_dcmp();
    void merge_mp();
    void merge_route(const std::vector<std::pair<int, int>>& a);
};

// Submodules
// trim - 최고차항의 계수가 0이 아닐 때까지 제거.
std::vector<double> trim_poly(
    const std::vector<double>& coeff
);

std::string vec_to_str(const std::vector<double>& v);

// attach - 두 복잡도의 complexity를 operation에 따라 결합
Complexity attach_poly(Complexity c1, Poly d2, Complexity c2, char attach_type);
Complexity mult_xi(XI d1, Poly d2, Complexity c2);

// 연산과정 클래스
struct step {
    char op;
    std::string key1;
    std::string key2;
    std::string save_key;
};

class EvalStep {
public:
    int term_count = 0;
    int coeff_count = 0;
    std::map<std::string, std::string> powers;
    std::map<std::string, std::string> coeffs;
    std::map<std::string, std::string> terms;
    std::deque<step> eval_step;

    EvalStep();
    EvalStep(std::shared_ptr<Decomp> dcmp);
    std::string serialize_dcmp(std::shared_ptr<Decomp> dcmp);
    std::string get_value(std::string key);
    std::string add_value(char key_header, int key_number, std::string data);
    std::string add_step(std::string key1, std::string key2, char op);
    void print_step();
    void print_line(step s, int& count);
};