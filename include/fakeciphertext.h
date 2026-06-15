#pragma once

#include <utility>
#include <vector>
#include "error_bound.h"

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
        Ciphertext(double x, int size, double err = 0.0, double scale = 1.0);
        Ciphertext(double a, double b, double interval_size = 1.0 / 16384.0, double err = 0.0, double scale = 1.0, bool encode = true);
        Ciphertext(const std::vector<std::pair<double, double>>& intervals, double err = 0.0, double scale = 1.0, bool encode = true);
        Ciphertext(std::vector<double>& x, double err = 0.0, double scale = 1.0, bool encode = true);
        Ciphertext(std::vector<double>& x, std::vector<double>& ct_high, std::vector<double>& ct_low, double err = 0.0, double scale = 1.0);

        void set_error(double err);
    };

    fct::Ciphertext add(fct::Ciphertext& ct1, fct::Ciphertext& ct2);
    fct::Ciphertext add_pt(fct::Ciphertext& ct, double scalar);
    fct::Ciphertext mult(fct::Ciphertext& ct1, fct::Ciphertext& ct2, ErrBound& eb);
    fct::Ciphertext square(fct::Ciphertext& ct, ErrBound& eb);
    fct::Ciphertext mult_pt(double scalar, fct::Ciphertext& ct, ErrBound& eb, bool rescale = true);

    std::pair<std::vector<double>, std::vector<double>> dec(fct::Ciphertext& ct);
    std::pair<double, double> dec_interval(fct::Ciphertext& ct);
    std::vector<std::pair<double, double>> dec_intervals(fct::Ciphertext& ct);

    fct::Ciphertext evaluate(fct::Ciphertext& x, EvalStep es, ErrBound& eb);
    fct::Ciphertext evaluate_composite(fct::Ciphertext& x, const std::vector<EvalStep>& eval_steps, ErrBound& eb);

    std::pair<double, double> estimate_range(ErrBound& eb, std::pair<double, double> interval, EvalStep& es, double interval_size, bool lipschitz = false);
    std::pair<double, double> estimate_composite_range(ErrBound& eb, std::pair<double, double> interval, const std::vector<EvalStep>& eval_steps, double interval_size, bool lipschitz = false);

    namespace Lipschitz
    {
        struct State
        {
            double pt_abs = 0.0;
            double high_abs = 0.0;
            double low_abs = 0.0;

            double pt_lip = 0.0;
            double high_lip = 0.0;
            double low_lip = 0.0;
        };

        struct Result
        {
            double M_L = 0.0;
            double M_U = 0.0;
            double M = 0.0;
        };

        Result estimate_lipschitz(ErrBound& eb, std::pair<double, double> interval, EvalStep& es, double input_err = 0.0);
        std::pair<double, double> estimate_ml_mu(ErrBound& eb, std::pair<double, double> interval, EvalStep& es, double input_err = 0.0);

        Result estimate_composite_lipschitz(ErrBound& eb, std::pair<double, double> interval, const std::vector<EvalStep>& eval_steps, double input_err = 0.0);
        std::pair<double, double> estimate_composite_ml_mu(ErrBound& eb, std::pair<double, double> interval, const std::vector<EvalStep>& eval_steps, double input_err = 0.0);

        /*
         * Composite Lipschitz estimation with stage-wise tightening.
         *
         * The old estimate_composite_lipschitz propagates only absolute-size
         * bounds through all composed polynomials. For indicator-function
         * routes, large alternating coefficients may cancel in the actual
         * range, but the absolute-size propagation cannot exploit this
         * cancellation and may produce unusably large M_L and M_U.
         *
         * This routine still propagates the Lipschitz constants through the
         * whole composite function, but after each polynomial stage it tightens
         * the magnitude bound using the ERE range of that stage. This keeps
         * the next stage's product-rule bound from being evaluated on an
         * artificial huge interval.
         */
        Result estimate_composite_lipschitz_tight(
            ErrBound& eb,
            std::pair<double, double> interval,
            const std::vector<EvalStep>& eval_steps,
            double interval_size,
            double input_err = 0.0
        );

        std::pair<double, double> estimate_composite_ml_mu_tight(
            ErrBound& eb,
            std::pair<double, double> interval,
            const std::vector<EvalStep>& eval_steps,
            double interval_size,
            double input_err = 0.0
        );
    }
}
