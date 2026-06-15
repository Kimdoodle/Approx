#pragma once
#include <vector>
#include "math.h"

// Plain data 전용
namespace pt
{
    std::vector<double> exp(const std::vector<double>& x, int exp);
    std::vector<double> negate(const std::vector<double>& x);
    std::vector<double> mult(const std::vector<double>& a, const std::vector<double>& b);
    std::vector<double> add(const std::vector<double>& a, const std::vector<double>& b);
    std::vector<double> sub(const std::vector<double>& a, const std::vector<double>& b);
    std::vector<double> mult_plain(const std::vector<double>& a, double scalar);
    std::vector<double> add_plain(const std::vector<double>& a, double scalar);
    std::vector<double> evaluate(const std::vector<double>& coeff, const std::vector<double>& x);
}