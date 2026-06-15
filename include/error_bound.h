#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <format>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "polyEval_class.h"
#include "eval_plain.h"


using INTERVAL = std::pair<std::pair<double, double>, std::pair<double, double>>;

class ErrBound 
{
public:
    double Bc, Bs, scale;

    ErrBound(double sigma, int N, int h, int s);

    // Bclean
    double B_clean(double sigma, int N, int h);

    // Bscale
    double B_scale(int N, int h);
};

std::string bit_diff(double a, double b);
std::string clearNum(double a);
