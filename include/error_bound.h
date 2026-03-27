#pragma once
#include <cmath>
#include <vector>
#include <utility>

#include "polyEval_class.h"
#include "eval_plain.h"
#include "CKKS_params.h"

class ErrBound 
{
public:
    double Bc, Bs, scale;
    ErrBound(double sigma, int N, int h, int s);
    
    // Bclean
    double B_clean(double sigma, int N, int h);

    // Bscale
    double B_scale(int N, int h);

    // Error bound of coeff.
    // double cal_bound(std::shared_ptr<Decomp> dcmp, std::vector<double>& x);
};

namespace fct
{
    class Ciphertext 
    {
    public:
        std::vector<double> pt;
        double err;
        std::vector<double> ct_high;
        std::vector<double> ct_low;
        double scale;

        Ciphertext();
        Ciphertext(double x, int size, double err=0.0, double scale=1.0);
        Ciphertext(std::vector<double>& x, double err=0.0, double scale=1.0, bool encode=true);
        Ciphertext(std::vector<double>& x, std::vector<double>& ct_high, std::vector<double>& ct_low, double err=0.0, double scale=1.0);
        void set_error(double err);
    };

    Ciphertext add(Ciphertext& ct1, Ciphertext& ct2);
    Ciphertext add_pt(Ciphertext& ct, double scalar);
    Ciphertext mult(Ciphertext& ct1, Ciphertext& ct2, ErrBound& eb);
    Ciphertext square(Ciphertext& ct, ErrBound& eb);
    Ciphertext mult_pt(double scalar, Ciphertext& ct, ErrBound& eb, bool rescale=true);

    std::pair<std::vector<double>, std::vector<double>> dec(Ciphertext& ct);
    Ciphertext evaluate(Ciphertext& x, EvalStep es, ErrBound& eb);
    Ciphertext evaluate_cl(Ciphertext& x, ErrBound& eb);
}

std::string bit_diff(double a, double b);
std::string clearNum(double a);