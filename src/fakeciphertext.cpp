#include "../include/fakeciphertext.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <format>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fct
{
    fct::Ciphertext::Ciphertext() = default;

    fct::Ciphertext::Ciphertext(double x, int size, double err, double scale)
    {
        std::vector<double> result_pt(size, x * scale);
        this->pt = result_pt;
        this->set_error(err);
        this->scale = scale;
    }

    fct::Ciphertext::Ciphertext(double a, double b, double interval_size, double err, double scale, bool encode)
    {
        int count = int(std::ceil((b - a) / interval_size)) + 1;

        this->pt.reserve(count);
        this->ct_high.reserve(count);
        this->ct_low.reserve(count);

        for(int i = 0; ; i++)
        {
            double mid = a + interval_size * double(i);

            if(mid > b)
                mid = b;

            if(encode)
                mid *= scale;

            this->pt.push_back(mid);
            this->ct_low.push_back(mid - err);
            this->ct_high.push_back(mid + err);

            if(mid == b || (encode && mid == b * scale))
                break;
        }

        this->err = err;
        this->scale = scale;
    }

    fct::Ciphertext::Ciphertext(const std::vector<std::pair<double, double>>& intervals, double err, double scale, bool encode)
    {
        this->pt.reserve(intervals.size());
        this->ct_high.reserve(intervals.size());
        this->ct_low.reserve(intervals.size());

        for(const auto& interval: intervals)
        {
            double low = interval.first;
            double high = interval.second;

            if(low > high)
                std::swap(low, high);

            double mid = (low + high) / 2.0;

            if(encode)
            {
                low *= scale;
                high *= scale;
                mid *= scale;
            }

            this->pt.push_back(mid);
            this->ct_low.push_back(low - err);
            this->ct_high.push_back(high + err);
        }

        this->err = err;
        this->scale = scale;
    }

    fct::Ciphertext::Ciphertext(std::vector<double>& x, double err, double scale, bool encode)
    {
        if(encode)
            this->pt = pt::mult_plain(x, scale);
        else
            this->pt = x;

        this->set_error(err);
        this->scale = scale;
    }

    fct::Ciphertext::Ciphertext(std::vector<double>& x, std::vector<double>& ct_high, std::vector<double>& ct_low, double err, double scale)
    {
        this->pt = x;
        this->err = err;
        this->ct_high = ct_high;
        this->ct_low = ct_low;
        this->scale = scale;
    }

    void fct::Ciphertext::set_error(double err)
    {
        this->err = err;
        this->ct_high = pt::add_plain(this->pt, err);
        this->ct_low = pt::add_plain(this->pt, -err);
    }

    fct::Ciphertext add(fct::Ciphertext& ct1, fct::Ciphertext& ct2)
    {
        assert(ct1.scale == ct2.scale);
        assert(ct1.pt.size() == ct2.pt.size());

        std::vector<double> result_pt = pt::add(ct1.pt, ct2.pt);
        std::vector<double> result_high = pt::add(ct1.ct_high, ct2.ct_high);
        std::vector<double> result_low = pt::add(ct1.ct_low, ct2.ct_low);

        return fct::Ciphertext(result_pt, result_high, result_low, 0.0, ct1.scale);
    }

    fct::Ciphertext add_pt(fct::Ciphertext& ct, double scalar)
    {
        scalar *= ct.scale;

        std::vector<double> result_pt = pt::add_plain(ct.pt, scalar);
        std::vector<double> result_high = pt::add_plain(ct.ct_high, scalar);
        std::vector<double> result_low = pt::add_plain(ct.ct_low, scalar);

        return fct::Ciphertext(result_pt, result_high, result_low, 0.0, ct.scale);
    }

    fct::Ciphertext mult(fct::Ciphertext& ct1, fct::Ciphertext& ct2, ErrBound& eb)
    {
        assert(ct1.pt.size() == ct2.pt.size());

        std::vector<double> result_pt = pt::mult(ct1.pt, ct2.pt);
        double err = eb.Bs;
        double scale = (ct1.scale * ct2.scale) / eb.scale;

        std::vector<double> bound1 = pt::mult(ct1.ct_high, ct2.ct_high);
        std::vector<double> bound2 = pt::mult(ct1.ct_high, ct2.ct_low);
        std::vector<double> bound3 = pt::mult(ct1.ct_low, ct2.ct_high);
        std::vector<double> bound4 = pt::mult(ct1.ct_low, ct2.ct_low);

        std::vector<double> result_high, result_low;
        result_high.reserve(result_pt.size());
        result_low.reserve(result_pt.size());

        for(size_t i = 0; i < result_pt.size(); i++)
        {
            result_high.push_back(std::max({bound1[i], bound2[i], bound3[i], bound4[i], result_pt[i]}));
            result_low.push_back(std::min({bound1[i], bound2[i], bound3[i], bound4[i], result_pt[i]}));
        }

        result_high = pt::mult_plain(result_high, 1.0 / eb.scale);
        result_low = pt::mult_plain(result_low, 1.0 / eb.scale);
        result_pt = pt::mult_plain(result_pt, 1.0 / eb.scale);

        result_high = pt::add_plain(result_high, err);
        result_low = pt::add_plain(result_low, -err);

        return fct::Ciphertext(result_pt, result_high, result_low, err, scale);
    }

    fct::Ciphertext square(fct::Ciphertext& ct, ErrBound& eb)
    {
        std::vector<double> result_pt = pt::mult(ct.pt, ct.pt);
        double err = eb.Bs;
        double scale = (ct.scale * ct.scale) / eb.scale;

        std::vector<double> result_high, result_low;
        result_high.reserve(result_pt.size());
        result_low.reserve(result_pt.size());

        for(size_t i = 0; i < result_pt.size(); i++)
        {
            const double low = ct.ct_low[i];
            const double high = ct.ct_high[i];
            const double low2 = low * low;
            const double high2 = high * high;

            double upper = std::max(low2, high2);
            double lower;

            if(low <= 0.0 && high >= 0.0)
                lower = 0.0;
            else
                lower = std::min(low2, high2);

            result_high.push_back(std::max({upper, result_pt[i]}));
            result_low.push_back(std::min({lower, result_pt[i]}));
        }

        result_high = pt::mult_plain(result_high, 1.0 / eb.scale);
        result_high = pt::add_plain(result_high, err);
        result_low = pt::mult_plain(result_low, 1.0 / eb.scale);
        result_low = pt::add_plain(result_low, -err);

        result_pt = pt::mult_plain(result_pt, 1.0 / eb.scale);

        return fct::Ciphertext(result_pt, result_high, result_low, err, scale);
    }

    fct::Ciphertext mult_pt(double scalar, fct::Ciphertext& ct, ErrBound& eb, bool rescale)
    {
        double err = 0.0;
        double scale = ct.scale;

        if(rescale)
        {
            scalar *= ct.scale;
            err = eb.Bs;
        }

        std::vector<double> result_pt = pt::mult_plain(ct.pt, scalar);
        std::vector<double> result_high = pt::mult_plain(ct.ct_high, scalar);
        std::vector<double> result_low = pt::mult_plain(ct.ct_low, scalar);

        if(scalar < 0.0)
            std::swap(result_high, result_low);

        if(rescale)
        {
            result_pt = pt::mult_plain(result_pt, 1.0 / eb.scale);
            result_high = pt::mult_plain(result_high, 1.0 / eb.scale);
            result_high = pt::add_plain(result_high, err);
            result_low = pt::mult_plain(result_low, 1.0 / eb.scale);
            result_low = pt::add_plain(result_low, -err);
        }

        return fct::Ciphertext(result_pt, result_high, result_low, err, scale);
    }

    std::pair<std::vector<double>, std::vector<double>> dec(fct::Ciphertext& ct)
    {
        std::vector<double> pt_high = pt::mult_plain(ct.ct_high, 1.0 / ct.scale);
        std::vector<double> pt_low = pt::mult_plain(ct.ct_low, 1.0 / ct.scale);

        return {pt_high, pt_low};
    }

    std::pair<double, double> dec_interval(fct::Ciphertext& ct)
    {
        auto [pt_high, pt_low] = dec(ct);

        double high = *std::max_element(pt_high.begin(), pt_high.end());
        double low = *std::min_element(pt_low.begin(), pt_low.end());

        return {low, high};
    }

    std::vector<std::pair<double, double>> dec_intervals(fct::Ciphertext& ct)
    {
        auto [pt_high, pt_low] = dec(ct);

        std::vector<std::pair<double, double>> intervals;
        intervals.reserve(pt_high.size());

        for(size_t i = 0; i < pt_high.size(); i++)
            intervals.push_back({pt_low[i], pt_high[i]});

        return intervals;
    }

    fct::Ciphertext evaluate(fct::Ciphertext& x, EvalStep es, ErrBound& eb)
    {
        fct::Ciphertext fctx = x;
        std::map<std::string, fct::Ciphertext> fpowers;
        std::map<std::string, double> coeffs;
        std::map<std::string, fct::Ciphertext> fterms;

        fpowers.insert_or_assign("P1", fctx);

        fct::Ciphertext res;

        for(auto step: es.eval_step)
        {
            switch(step.op)
            {
                case 'o':
                    if(step.key1[0] != 'C')
                    {
                        if(step.key1[0] == 'P')
                            continue;

                        throw std::invalid_argument(std::format("Main::RES: invalid key data {} {}", step.key1, step.key2));
                    }

                    coeffs.insert_or_assign(step.key1, std::stod(step.key2));
                    break;

                case '+':
                    if(step.key1[0] == 'T' && step.key2[0] == 'C')
                    {
                        res = add_pt(fterms.at(step.key1), coeffs.at(step.key2));
                    }
                    else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                    {
                        res = add(fterms.at(step.key1), fterms.at(step.key2));
                    }
                    else if(step.key1[0] == 'P' && step.key2[0] == 'C')
                    {
                        res = add_pt(fpowers.at(step.key1), coeffs.at(step.key2));
                    }
                    else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                    {
                        res = add(fpowers.at(step.key1), fterms.at(step.key2));
                    }
                    else if(step.key1[0] == 'T' && step.key2[0] == 'P')
                    {
                        res = add(fterms.at(step.key1), fpowers.at(step.key2));
                    }
                    else
                    {
                        throw std::invalid_argument(std::format("Main::ADD: invalid key data {} {}", step.key1, step.key2));
                    }

                    fterms.insert_or_assign(step.save_key, res);
                    break;

                case 'x':
                    if(step.key1[0] == 'P' && step.key2[0] == 'P')
                    {
                        if(step.key1 == step.key2)
                            res = square(fpowers.at(step.key1), eb);
                        else
                            res = mult(fpowers.at(step.key1), fpowers.at(step.key2), eb);
                    }
                    else if(step.key1[0] == 'C' && step.key2[0] == 'P')
                    {
                        res = mult_pt(coeffs.at(step.key1), fpowers.at(step.key2), eb);
                    }
                    else if(step.key1[0] == 'P' && step.key2[0] == 'C')
                    {
                        res = mult_pt(coeffs.at(step.key2), fpowers.at(step.key1), eb, true);
                    }
                    else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                    {
                        res = mult(fpowers.at(step.key1), fterms.at(step.key2), eb);
                    }
                    else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                    {
                        res = mult(fterms.at(step.key1), fterms.at(step.key2), eb);
                    }
                    else if(step.key1[0] == 'T' && step.key2[0] == 'P')
                    {
                        res = mult(fterms.at(step.key1), fpowers.at(step.key2), eb);
                    }
                    else
                    {
                        throw std::invalid_argument(std::format("Main::MUL: invalid key data {} {}", step.key1, step.key2));
                    }

                    switch(step.save_key[0])
                    {
                        case 'P':
                            fpowers.insert_or_assign(step.save_key, res);
                            break;

                        case 'T':
                            fterms.insert_or_assign(step.save_key, res);
                            break;

                        default:
                            break;
                    }
                    break;

                default:
                    break;
            }
        }

        fct::Ciphertext fres = fterms.at(es.eval_step.back().save_key);
        return fres;
    }

    fct::Ciphertext evaluate_composite(fct::Ciphertext& x, const std::vector<EvalStep>& eval_steps, ErrBound& eb)
    {
        if(eval_steps.empty())
            throw std::invalid_argument("evaluate_composite: empty evaluation step list");

        fct::Ciphertext res = x;

        for(const auto& es: eval_steps)
            res = evaluate(res, es, eb);

        return res;
    }

    std::pair<double, double> estimate_range(ErrBound& eb, std::pair<double, double> interval, EvalStep& es, double interval_size, bool lipschitz)
    {
        std::pair<double, double> range;

        // const double interval_size = 1.0 / 16384.0;

        fct::Ciphertext ct = fct::Ciphertext(interval.first, interval.second, interval_size, 0.0, eb.scale, true);

        ct = evaluate(ct, es, eb);

        // debug - 구간 크기 측정
        std::vector<double> slot_size(ct.pt.size());

        for(size_t i = 0; i < ct.pt.size(); i++)
        {
            double low = ct.ct_low[i] / ct.scale;
            double high = ct.ct_high[i] / ct.scale;
            slot_size[i] = high - low;
        }

        double sum = std::accumulate(slot_size.begin(), slot_size.end(), 0.0);
        // std::cout << "Slot size: " << clearNum(sum / static_cast<double>(slot_size.size())) << std::endl;

        range = dec_interval(ct);

        if(lipschitz)
        {
            auto lip = Lipschitz::estimate_lipschitz(eb, interval, es);

            const double lower_padding = 0.5 * interval_size * lip.M_L;
            const double upper_padding = 0.5 * interval_size * lip.M_U;

            range.first -= lower_padding;
            range.second += upper_padding;

            // std::cout << "Lipschitz M_L: " << clearNum(lip.M_L) << std::endl;
            // std::cout << "Lipschitz M_U: " << clearNum(lip.M_U) << std::endl;
            // std::cout << "Lipschitz lower padding: "
            //         << clearNum(lower_padding) << std::endl;
            // std::cout << "Lipschitz upper padding: "
            //         << clearNum(upper_padding) << std::endl;
        }
        return range;
    }

    namespace Lipschitz
    {
        namespace
        {
            double max_abs_interval(double a, double b)
            {
                if(a > b)
                    std::swap(a, b);

                return std::max(std::abs(a), std::abs(b));
            }

            double prod_lip(double a_abs, double a_lip, double b_abs, double b_lip)
            {
                return a_lip * b_abs + a_abs * b_lip;
            }

            State make_input(std::pair<double, double> interval, double input_err, double scale)
            {
                double a = interval.first;
                double b = interval.second;

                if(a > b)
                    std::swap(a, b);

                const double beta_clean = input_err / scale;

                State res;
                res.pt_abs = max_abs_interval(a, b);
                res.high_abs = max_abs_interval(a + beta_clean, b + beta_clean);
                res.low_abs = max_abs_interval(a - beta_clean, b - beta_clean);

                res.pt_lip = 1.0;
                res.high_lip = 1.0;
                res.low_lip = 1.0;

                return res;
            }

            State make_const(double scalar)
            {
                State res;
                res.pt_abs = std::abs(scalar);
                res.high_abs = std::abs(scalar);
                res.low_abs = std::abs(scalar);
                return res;
            }

            State add_state(const State& lhs, const State& rhs)
            {
                State res;
                res.pt_abs = lhs.pt_abs + rhs.pt_abs;
                res.high_abs = lhs.high_abs + rhs.high_abs;
                res.low_abs = lhs.low_abs + rhs.low_abs;

                res.pt_lip = lhs.pt_lip + rhs.pt_lip;
                res.high_lip = lhs.high_lip + rhs.high_lip;
                res.low_lip = lhs.low_lip + rhs.low_lip;

                return res;
            }

            State add_plain_state(const State& ct, double scalar)
            {
                State res = ct;

                // |g(x) + c| <= |g(x)| + |c|.  The Lipschitz constant is unchanged.
                const double abs_scalar = std::abs(scalar);
                res.pt_abs += abs_scalar;
                res.high_abs += abs_scalar;
                res.low_abs += abs_scalar;

                return res;
            }

            State mult_plain_state(double scalar, const State& ct, ErrBound& eb, bool rescale)
            {
                State res;
                const double abs_scalar = std::abs(scalar);
                const double beta_scale = rescale ? eb.Bs / eb.scale : 0.0;

                res.pt_abs = abs_scalar * ct.pt_abs;
                res.pt_lip = abs_scalar * ct.pt_lip;

                if(scalar >= 0.0)
                {
                    res.high_abs = abs_scalar * ct.high_abs + beta_scale;
                    res.low_abs = abs_scalar * ct.low_abs + beta_scale;
                    res.high_lip = abs_scalar * ct.high_lip;
                    res.low_lip = abs_scalar * ct.low_lip;
                }
                else
                {
                    // Multiplication by a negative scalar swaps the upper/lower envelopes.
                    res.high_abs = abs_scalar * ct.low_abs + beta_scale;
                    res.low_abs = abs_scalar * ct.high_abs + beta_scale;
                    res.high_lip = abs_scalar * ct.low_lip;
                    res.low_lip = abs_scalar * ct.high_lip;
                }

                return res;
            }

            State square_state(const State& ct, ErrBound& eb)
            {
                State res;
                const double beta_scale = eb.Bs / eb.scale;

                const double pt_abs = ct.pt_abs * ct.pt_abs;
                const double high_abs = ct.high_abs * ct.high_abs;
                const double low_abs = ct.low_abs * ct.low_abs;

                const double pt_lip = 2.0 * ct.pt_abs * ct.pt_lip;
                const double high_lip = 2.0 * ct.high_abs * ct.high_lip;
                const double low_lip = 2.0 * ct.low_abs * ct.low_lip;

                const double envelope_abs = std::max({pt_abs, high_abs, low_abs});
                const double envelope_lip = std::max({pt_lip, high_lip, low_lip, 0.0});

                res.pt_abs = pt_abs;
                res.pt_lip = pt_lip;

                // upper = max(low^2, high^2, pt^2) + beta_scale
                // lower = min(low^2, high^2, pt^2, 0 if the interval crosses zero) - beta_scale
                // min/max of M-Lipschitz functions is M-Lipschitz with M equal to the max of constants.
                res.high_abs = envelope_abs + beta_scale;
                res.low_abs = envelope_abs + beta_scale;
                res.high_lip = envelope_lip;
                res.low_lip = envelope_lip;

                return res;
            }

            State mult_state(const State& lhs, const State& rhs, ErrBound& eb)
            {
                State res;
                const double beta_scale = eb.Bs / eb.scale;

                res.pt_abs = lhs.pt_abs * rhs.pt_abs;
                res.pt_lip = prod_lip(lhs.pt_abs, lhs.pt_lip, rhs.pt_abs, rhs.pt_lip);

                const double hh_abs = lhs.high_abs * rhs.high_abs;
                const double hl_abs = lhs.high_abs * rhs.low_abs;
                const double lh_abs = lhs.low_abs * rhs.high_abs;
                const double ll_abs = lhs.low_abs * rhs.low_abs;

                const double hh_lip = prod_lip(lhs.high_abs, lhs.high_lip, rhs.high_abs, rhs.high_lip);
                const double hl_lip = prod_lip(lhs.high_abs, lhs.high_lip, rhs.low_abs, rhs.low_lip);
                const double lh_lip = prod_lip(lhs.low_abs, lhs.low_lip, rhs.high_abs, rhs.high_lip);
                const double ll_lip = prod_lip(lhs.low_abs, lhs.low_lip, rhs.low_abs, rhs.low_lip);

                const double envelope_abs = std::max({res.pt_abs, hh_abs, hl_abs, lh_abs, ll_abs});
                const double envelope_lip = std::max({res.pt_lip, hh_lip, hl_lip, lh_lip, ll_lip});

                // upper/lower are max/min over endpoint products and the nominal product, followed by rescale padding.
                res.high_abs = envelope_abs + beta_scale;
                res.low_abs = envelope_abs + beta_scale;
                res.high_lip = envelope_lip;
                res.low_lip = envelope_lip;

                return res;
            }

            const State& get_state(const std::string& key, const std::map<std::string, State>& fpowers, const std::map<std::string, State>& fterms)
            {
                if(key.empty())
                    throw std::invalid_argument("Lipschitz::get_state: empty key");

                if(key[0] == 'P')
                    return fpowers.at(key);

                if(key[0] == 'T')
                    return fterms.at(key);

                throw std::invalid_argument(std::format("Lipschitz::get_state: invalid key {}", key));
            }

            State evaluate_state(const State& input, ErrBound& eb, const EvalStep& es)
            {
                std::map<std::string, State> fpowers;
                std::map<std::string, double> coeffs;
                std::map<std::string, State> fterms;

                fpowers.insert_or_assign("P1", input);

                State res;

                for(const auto& step: es.eval_step)
                {
                    switch(step.op)
                    {
                        case 'o':
                            if(step.key1[0] != 'C')
                            {
                                if(step.key1[0] == 'P')
                                    continue;

                                throw std::invalid_argument(std::format("Lipschitz::RES: invalid key data {} {}", step.key1, step.key2));
                            }

                            coeffs.insert_or_assign(step.key1, std::stod(step.key2));
                            break;

                        case '+':
                            if(step.key1[0] == 'T' && step.key2[0] == 'C')
                            {
                                res = add_plain_state(fterms.at(step.key1), coeffs.at(step.key2));
                            }
                            else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                            {
                                res = add_state(fterms.at(step.key1), fterms.at(step.key2));
                            }
                            else if(step.key1[0] == 'P' && step.key2[0] == 'C')
                            {
                                res = add_plain_state(fpowers.at(step.key1), coeffs.at(step.key2));
                            }
                            else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                            {
                                res = add_state(fpowers.at(step.key1), fterms.at(step.key2));
                            }
                            else if(step.key1[0] == 'T' && step.key2[0] == 'P')
                            {
                                res = add_state(fterms.at(step.key1), fpowers.at(step.key2));
                            }
                            else
                            {
                                throw std::invalid_argument(std::format("Lipschitz::ADD: invalid key data {} {}", step.key1, step.key2));
                            }

                            fterms.insert_or_assign(step.save_key, res);
                            break;

                        case 'x':
                            if(step.key1[0] == 'P' && step.key2[0] == 'P')
                            {
                                if(step.key1 == step.key2)
                                    res = square_state(fpowers.at(step.key1), eb);
                                else
                                    res = mult_state(fpowers.at(step.key1), fpowers.at(step.key2), eb);
                            }
                            else if(step.key1[0] == 'C' && step.key2[0] == 'P')
                            {
                                res = mult_plain_state(coeffs.at(step.key1), fpowers.at(step.key2), eb, true);
                            }
                            else if(step.key1[0] == 'P' && step.key2[0] == 'C')
                            {
                                res = mult_plain_state(coeffs.at(step.key2), fpowers.at(step.key1), eb, true);
                            }
                            else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                            {
                                res = mult_state(fpowers.at(step.key1), fterms.at(step.key2), eb);
                            }
                            else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                            {
                                res = mult_state(fterms.at(step.key1), fterms.at(step.key2), eb);
                            }
                            else if(step.key1[0] == 'T' && step.key2[0] == 'P')
                            {
                                res = mult_state(fterms.at(step.key1), fpowers.at(step.key2), eb);
                            }
                            else
                            {
                                throw std::invalid_argument(std::format("Lipschitz::MUL: invalid key data {} {}", step.key1, step.key2));
                            }

                            switch(step.save_key[0])
                            {
                                case 'P':
                                    fpowers.insert_or_assign(step.save_key, res);
                                    break;

                                case 'T':
                                    fterms.insert_or_assign(step.save_key, res);
                                    break;

                                default:
                                    break;
                            }
                            break;

                        default:
                            break;
                    }
                }

                if(es.eval_step.empty())
                    throw std::invalid_argument("Lipschitz::evaluate_state: empty evaluation step");

                const std::string& output_key = es.eval_step.back().save_key;
                return get_state(output_key, fpowers, fterms);
            }

            State tighten_abs_with_interval(const State& state, std::pair<double, double> interval)
            {
                State res = state;

                const double abs_bound = max_abs_interval(interval.first, interval.second);

                /*
                 * Preserve the derivative/Lipschitz information, but replace
                 * the magnitude bound by the actual ERE range obtained after
                 * the current polynomial stage.
                 *
                 * This is the crucial step. Without this reset, an alternating
                 * polynomial such as a0 - a2 x^2 + a4 x^4 - ... is treated as
                 * if all absolute terms could add constructively. The next
                 * polynomial then receives an artificial input magnitude and
                 * the product-rule bound explodes.
                 */
                res.pt_abs = abs_bound;
                res.high_abs = abs_bound;
                res.low_abs = abs_bound;

                return res;
            }

            Result make_result(const State& output)
            {
                Result result;
                result.M_L = output.low_lip;
                result.M_U = output.high_lip;
                result.M = std::max(result.M_L, result.M_U);
                return result;
            }
        }

        Result estimate_lipschitz(ErrBound& eb, std::pair<double, double> interval, EvalStep& es, double input_err)
        {
            std::map<std::string, State> fpowers;
            std::map<std::string, double> coeffs;
            std::map<std::string, State> fterms;

            fpowers.insert_or_assign("P1", make_input(interval, input_err, eb.scale));

            State res;

            for(const auto& step: es.eval_step)
            {
                switch(step.op)
                {
                    case 'o':
                        if(step.key1[0] != 'C')
                        {
                            if(step.key1[0] == 'P')
                                continue;

                            throw std::invalid_argument(std::format("Lipschitz::RES: invalid key data {} {}", step.key1, step.key2));
                        }

                        coeffs.insert_or_assign(step.key1, std::stod(step.key2));
                        break;

                    case '+':
                        if(step.key1[0] == 'T' && step.key2[0] == 'C')
                        {
                            res = add_plain_state(fterms.at(step.key1), coeffs.at(step.key2));
                        }
                        else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                        {
                            res = add_state(fterms.at(step.key1), fterms.at(step.key2));
                        }
                        else if(step.key1[0] == 'P' && step.key2[0] == 'C')
                        {
                            res = add_plain_state(fpowers.at(step.key1), coeffs.at(step.key2));
                        }
                        else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                        {
                            res = add_state(fpowers.at(step.key1), fterms.at(step.key2));
                        }
                        else if(step.key1[0] == 'T' && step.key2[0] == 'P')
                        {
                            res = add_state(fterms.at(step.key1), fpowers.at(step.key2));
                        }
                        else
                        {
                            throw std::invalid_argument(std::format("Lipschitz::ADD: invalid key data {} {}", step.key1, step.key2));
                        }

                        fterms.insert_or_assign(step.save_key, res);
                        break;

                    case 'x':
                        if(step.key1[0] == 'P' && step.key2[0] == 'P')
                        {
                            if(step.key1 == step.key2)
                                res = square_state(fpowers.at(step.key1), eb);
                            else
                                res = mult_state(fpowers.at(step.key1), fpowers.at(step.key2), eb);
                        }
                        else if(step.key1[0] == 'C' && step.key2[0] == 'P')
                        {
                            res = mult_plain_state(coeffs.at(step.key1), fpowers.at(step.key2), eb, true);
                        }
                        else if(step.key1[0] == 'P' && step.key2[0] == 'C')
                        {
                            res = mult_plain_state(coeffs.at(step.key2), fpowers.at(step.key1), eb, true);
                        }
                        else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                        {
                            res = mult_state(fpowers.at(step.key1), fterms.at(step.key2), eb);
                        }
                        else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                        {
                            res = mult_state(fterms.at(step.key1), fterms.at(step.key2), eb);
                        }
                        else if(step.key1[0] == 'T' && step.key2[0] == 'P')
                        {
                            res = mult_state(fterms.at(step.key1), fpowers.at(step.key2), eb);
                        }
                        else
                        {
                            throw std::invalid_argument(std::format("Lipschitz::MUL: invalid key data {} {}", step.key1, step.key2));
                        }

                        switch(step.save_key[0])
                        {
                            case 'P':
                                fpowers.insert_or_assign(step.save_key, res);
                                break;

                            case 'T':
                                fterms.insert_or_assign(step.save_key, res);
                                break;

                            default:
                                break;
                        }
                        break;

                    default:
                        break;
                }
            }

            if(es.eval_step.empty())
                throw std::invalid_argument("Lipschitz::estimate_lipschitz: empty evaluation step");

            const std::string& output_key = es.eval_step.back().save_key;
            const State& output = get_state(output_key, fpowers, fterms);

            Result result;
            result.M_L = output.low_lip;
            result.M_U = output.high_lip;
            result.M = std::max(result.M_L, result.M_U);

            return result;
        }

        std::pair<double, double> estimate_ml_mu(ErrBound& eb, std::pair<double, double> interval, EvalStep& es, double input_err)
        {
            Result result = estimate_lipschitz(eb, interval, es, input_err);
            return {result.M_L, result.M_U};
        }

        Result estimate_composite_lipschitz(ErrBound& eb, std::pair<double, double> interval, const std::vector<EvalStep>& eval_steps, double input_err)
        {
            if(eval_steps.empty())
                throw std::invalid_argument("Lipschitz::estimate_composite_lipschitz: empty evaluation step list");

            State current = make_input(interval, input_err, eb.scale);

            for(const auto& es: eval_steps)
                current = evaluate_state(current, eb, es);

            return make_result(current);
        }

        std::pair<double, double> estimate_composite_ml_mu(ErrBound& eb, std::pair<double, double> interval, const std::vector<EvalStep>& eval_steps, double input_err)
        {
            Result result = estimate_composite_lipschitz(eb, interval, eval_steps, input_err);
            return {result.M_L, result.M_U};
        }

        Result estimate_composite_lipschitz_tight(
            ErrBound& eb,
            std::pair<double, double> interval,
            const std::vector<EvalStep>& eval_steps,
            double interval_size,
            double input_err
        )
        {
            if(eval_steps.empty())
                throw std::invalid_argument("Lipschitz::estimate_composite_lipschitz_tight: empty evaluation step list");

            if(interval_size <= 0.0)
                throw std::invalid_argument("Lipschitz::estimate_composite_lipschitz_tight: interval_size must be positive");

            State current = make_input(interval, input_err, eb.scale);
            std::pair<double, double> current_interval = interval;

            for(const auto& es: eval_steps)
            {
                /*
                 * 1. Propagate derivative bounds through the current stage.
                 *    The Lipschitz constants in `stage_state` are still with
                 *    respect to the original input because `current` already
                 *    carries the accumulated Lipschitz constants.
                 */
                State stage_state = evaluate_state(current, eb, es);

                /*
                 * 2. Compute a tight ERE range for the same stage.
                 *    This is an internal C++ call; Python does not need to
                 *    implement any ERE logic.
                 */
                EvalStep es_copy = es;
                std::pair<double, double> sampled_range =
                    fct::estimate_range(eb, current_interval, es_copy, interval_size, false);

                /*
                 * 3. Certify the stage range by the current stage Lipschitz
                 *    margin. This certified interval is used only to tighten
                 *    the magnitude bound for later stages.
                 */
                const double lower_padding = 0.5 * interval_size * stage_state.low_lip;
                const double upper_padding = 0.5 * interval_size * stage_state.high_lip;

                std::pair<double, double> certified_stage_range = {
                    sampled_range.first - lower_padding,
                    sampled_range.second + upper_padding
                };

                /*
                 * 4. Reset only the magnitude bound. Do not reset the
                 *    Lipschitz constants, because they represent the
                 *    accumulated sensitivity of the current composite prefix.
                 */
                current = tighten_abs_with_interval(stage_state, certified_stage_range);
                current_interval = certified_stage_range;
            }

            return make_result(current);
        }

        std::pair<double, double> estimate_composite_ml_mu_tight(
            ErrBound& eb,
            std::pair<double, double> interval,
            const std::vector<EvalStep>& eval_steps,
            double interval_size,
            double input_err
        )
        {
            Result result = estimate_composite_lipschitz_tight(eb, interval, eval_steps, interval_size, input_err);
            return {result.M_L, result.M_U};
        }
    }

    std::pair<double, double> estimate_composite_range(ErrBound& eb, std::pair<double, double> interval, const std::vector<EvalStep>& eval_steps, double interval_size, bool lipschitz)
    {
        if(eval_steps.empty())
            throw std::invalid_argument("estimate_composite_range: empty evaluation step list");

        fct::Ciphertext ct = fct::Ciphertext(interval.first, interval.second, interval_size, 0.0, eb.scale, true);
        ct = evaluate_composite(ct, eval_steps, eb);

        std::pair<double, double> range = dec_interval(ct);

        if(lipschitz)
        {
            auto lip = Lipschitz::estimate_composite_lipschitz_tight(eb, interval, eval_steps, interval_size);

            const double lower_padding = 0.5 * interval_size * lip.M_L;
            const double upper_padding = 0.5 * interval_size * lip.M_U;

            range.first -= lower_padding;
            range.second += upper_padding;
        }

        return range;
    }
}
