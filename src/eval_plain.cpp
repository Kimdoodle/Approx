#include "../include/eval_plain.h"

// Operations of Plain vector.
std::vector<double> pt::exp(const std::vector<double>& x, int exp)
{
    std::vector<double> res(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        res[i] = pow(x[i], exp);
    }
    return res;
}

std::vector<double> pt::negate(const std::vector<double>& x)
{
    return mult_plain(x, -1.0);
}

std::vector<double> pt::mult(const std::vector<double>& a, const std::vector<double>& b)
{
    std::vector<double> res(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        res[i] = a[i] * b[i];
    }
    return res;
}

std::vector<double> pt::add(const std::vector<double>& a, const std::vector<double>& b)
{
    std::vector<double> res(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        res[i] = a[i] + b[i];
    }
    return res;
}

std::vector<double> pt::sub(const std::vector<double>& a, const std::vector<double>& b)
{
    std::vector<double> res(a.size());
    res = add(a, negate(b));
    return res;
}

std::vector<double> pt::mult_plain(const std::vector<double>& a, double scalar)
{
    std::vector<double> res(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        res[i] = a[i] * scalar;
    }
    return res;
}

std::vector<double> pt::add_plain(const std::vector<double>& a, double scalar)
{
    std::vector<double> res(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        res[i] = a[i] + scalar;
    }
    return res;
}

std::vector<double> pt::eval(const std::vector<double>& coeff, const std::vector<double>& x)
{
    std::vector<double> res(x.size(), 0.0);

    if (coeff.empty()) return res;

    for (std::size_t i = 0; i < x.size(); ++i) {
        double acc = 0.0;
        for (std::size_t j = coeff.size(); j-- > 0; ) {
            acc = acc * x[i] + coeff[j];
        }
        res[i] = acc;
    }
    return res;
}
