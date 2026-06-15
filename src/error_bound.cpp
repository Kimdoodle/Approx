#include "../include/error_bound.h"

ErrBound::ErrBound(double sigma, int N, int h, int s)
{
    Bc = B_clean(sigma, N, h);
    Bs = B_scale(N, h);
    scale = pow(2, s);
}

double ErrBound::B_clean(double sigma, int N, int h) 
{
    return (8 * sqrt(2.0) * sigma * N)
            + (6 * sigma * sqrt(static_cast<double>(N)))
            + (16 * sigma * sqrt(static_cast<double>(h) * N));
}

double ErrBound::B_scale(int N, int h) 
{
    return sqrt(N / 3.0) * (3.0 + 8.0 * sqrt(static_cast<double>(h)));
}

std::string bit_diff(double bound, double data) 
{
    double diff = std::abs(bound - data);
    char sign = bound > data ? '-' : '+';

    if(diff == 0.0)
        return "0 bits";

    double log2_diff = std::abs(std::log2(diff));

    return std::format("{}{:.4f} bits", sign, log2_diff);
}

std::string clearNum(double a) 
{
    int b = static_cast<int>(std::round(a));
    char sign = a > b ? '+' : '-';

    const double e = std::abs(a - static_cast<double>(b));

    if(e == 0.0)
        return std::to_string(b);

    const double log2e = std::log2(e);

    return std::format("{} {}2^({:.10f})", b, sign, log2e);
}
