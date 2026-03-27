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
    // return sqrt(N * 3.0) + (8 * sqrt((static_cast<double>(h) * N) / 3.0));
}


/* double ErrBound::cal_bound(std::shared_ptr<Decomp> dcmp, std::vector<double>& x)
{    
    // 평문상태의 연산수행 결과
    std::vector<double> result_pt = pt::eval(dcmp->coeff, x);

    // 암호문상태의 연산결과 추정
    std::vector<double> result_ct;
    result_ct.reserve(result_pt.size());
    
    // 실제 연산과 동일한 과정을 수행함.
    // 1. x^i들 연산
    std::vector<double> ctx = x;
    fct::Ciphertext pt = fct::Ciphertext(ctx, scale);
    pt.set_error(0.0);
    std::map<int, fct::Ciphertext> powers = {{1, pt}};

    for(std::pair<int, int> route: dcmp->total_routes)
    {
        int a = route.first;
        int b = route.second;
        if(!dcmp->made_powers.contains(a+b))
            continue;
        fct::Ciphertext cta = powers.at(a);
        fct::Ciphertext ctb = powers.at(b);
        fct::Ciphertext ctres = fct::mult(cta, ctb, *this);
        powers.insert({a+b, ctres});
    }
    return 0.0;
} */

namespace fct
{
    Ciphertext::Ciphertext() = default;
    Ciphertext::Ciphertext(double x, int size, double err, double scale)
    {
        std::vector<double> pt(size, x*scale);
        this->pt = pt;
        this->set_error(err);
        this->scale = scale;
    }
    Ciphertext::Ciphertext(std::vector<double>& x, double err, double scale, bool encode)
    {
        if(encode)
            this->pt = pt::mult_plain(x, scale);
        else
            this->pt = x;
        this->set_error(err);
        this->scale = scale;
    }
    Ciphertext::Ciphertext(std::vector<double>& x, std::vector<double>& ct_high, std::vector<double>& ct_low, double err, double scale)
    {
        this->pt = x;
        this->err = err;
        this->ct_high = ct_high;
        this->ct_low = ct_low;
        this->scale = scale;
    }

    void Ciphertext::set_error(double err)
    {
        this->err = err;
        this->ct_high = pt::add_plain(this->pt, err);
        this->ct_low = pt::add_plain(this->pt, -err);
    }

    Ciphertext add(Ciphertext& ct1, Ciphertext& ct2)
    {
        assert(ct1.scale == ct2.scale);

        std::vector<double> pt = pt::add(ct1.pt, ct2.pt);
        std::vector<double> ct_high = pt::add(ct1.ct_high, ct2.ct_high);
        std::vector<double> ct_low = pt::add(ct1.ct_low, ct2.ct_low);
        return Ciphertext(pt, ct_high, ct_low, 0.0, ct1.scale);
    }

    Ciphertext add_pt(Ciphertext& ct, double scalar)
    {
        scalar *= ct.scale;
        std::vector<double> pt = pt::add_plain(ct.pt, scalar);
        std::vector<double> ct_high = pt::add_plain(ct.ct_high, scalar);
        std::vector<double> ct_low = pt::add_plain(ct.ct_low, scalar);
        return Ciphertext(pt, ct_high, ct_low, 0.0, ct.scale);
    }

    Ciphertext mult(Ciphertext& ct1, Ciphertext& ct2, ErrBound& eb)
    {
        std::vector<double> pt = pt::mult(ct1.pt, ct2.pt);
        double err = eb.Bs;
        double scale = (ct1.scale * ct2.scale) / eb.scale;

        // ct_high, ct_low 연산
        // 두 ct의 ct_high, ct_low의 hardmard product를 계산한 후 최댓값/최솟값을 선정.
        std::vector<double> bound1 = pt::mult(ct1.ct_high, ct2.ct_high);
        // std::vector<double> bound2 = mult_vector(ct1.ct_high, ct2.ct_low);
        // std::vector<double> bound3 = mult_vector(ct1.ct_low, ct2.ct_high);
        std::vector<double> bound4 = pt::mult(ct1.ct_low, ct2.ct_low);

        std::vector<double> ct_high, ct_low;
        ct_high.reserve(pt.size());
        ct_low.reserve(pt.size());
        
        for(int i=0; i<pt.size(); i++)
        {
            ct_high.push_back(std::max({bound1[i], bound4[i]}));
            ct_low.push_back(std::min({bound1[i], bound4[i]}));
        }

        ct_high = pt::mult_plain(ct_high, 1/eb.scale);
        ct_low = pt::mult_plain(ct_low, 1/eb.scale);
        pt = pt::mult_plain(pt, 1/eb.scale);
        
        ct_high = pt::add_plain(ct_high, err);
        ct_low = pt::add_plain(ct_low, -err);

        return Ciphertext(pt, ct_high, ct_low, err, scale);
    }

    Ciphertext square(Ciphertext& ct, ErrBound& eb)
    {
        std::vector<double> pt = pt::mult(ct.pt, ct.pt);
        double err = eb.Bs;
        double scale = (ct.scale * ct.scale) / eb.scale;

        // ct_high, ct_low 연산
        std::vector<double> bound1 = pt::mult(ct.ct_high, ct.ct_high);
        std::vector<double> bound4 = pt::mult(ct.ct_low, ct.ct_low);

        std::vector<double> ct_high, ct_low;
        ct_high.reserve(pt.size());
        ct_low.reserve(pt.size());
        
        for(int i=0; i<pt.size(); i++)
        {
            ct_high.push_back(std::max({bound1[i], pt[i], bound4[i]}));
            ct_low.push_back(std::min({bound1[i], pt[i], bound4[i]}));
        }

        ct_high = pt::mult_plain(ct_high, 1/eb.scale);
        ct_high = pt::add_plain(ct_high, err);
        ct_low = pt::mult_plain(ct_low, 1/eb.scale);
        ct_low = pt::add_plain(ct_low, -err);

        pt = pt::mult_plain(pt, 1/eb.scale);

        return Ciphertext(pt, ct_high, ct_low, err, scale);
    }

    Ciphertext mult_pt(double scalar, Ciphertext&ct, ErrBound& eb, bool rescale)
    {
        double err = 0.0;
        double scale = ct.scale;
        if(rescale)
        {
            scalar *= ct.scale;
            err = eb.Bs;
        }

        std::vector<double> pt = pt::mult_plain(ct.pt, scalar);
        std::vector<double> ct_high = pt::mult_plain(ct.ct_high, scalar);
        std::vector<double> ct_low = pt::mult_plain(ct.ct_low, scalar);


        if(scalar < 0.0)
        {
            std::swap(ct_high, ct_low);
        }
        if(rescale)
        {
            pt = pt::mult_plain(pt, 1.0/eb.scale);
            ct_high = pt::mult_plain(ct_high, 1.0/eb.scale);
            ct_high = pt::add_plain(ct_high, err);
            ct_low = pt::mult_plain(ct_low, 1.0/eb.scale);
            ct_low = pt::add_plain(ct_low, -err);
        }

        return Ciphertext(pt, ct_high, ct_low, err, scale);
    }

    std::pair<std::vector<double>, std::vector<double>> dec(Ciphertext& ct)
    {
        std::vector<double> pt_high = pt::mult_plain(ct.ct_high, 1.0/ct.scale);
        std::vector<double> pt_low = pt::mult_plain(ct.ct_low, 1.0/ct.scale);
        return {pt_high, pt_low};
    }

    Ciphertext evaluate(Ciphertext& x, EvalStep es, ErrBound& eb)
    {
        fct::Ciphertext fctx = x;
        std::map<std::string, fct::Ciphertext> fpowers;
        std::map<std::string, double> coeffs;
        std::map<std::string, fct::Ciphertext> fterms;
        fpowers.insert_or_assign("P1", fctx);
        fct::Ciphertext fvalue1, fvalue2, res;

        // 평가
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
                        res = fct::add_pt(fterms.at(step.key1), coeffs.at(step.key2));
                    }
                    else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                    {
                        res = fct::add(fterms.at(step.key1), fterms.at(step.key2));
                    }
                    else if(step.key1[0] == 'P' && step.key2[0] == 'C')
                    {
                        res = fct::add_pt(fpowers.at(step.key1), coeffs.at(step.key2));
                    }
                    else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                    {
                        res = fct::add(fpowers.at(step.key1), fterms.at(step.key2));
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
                        res = fct::mult(fpowers.at(step.key1), fpowers.at(step.key2), eb);
                    }
                    else if(step.key1[0] == 'C' && step.key2[0] == 'P')
                    {
                        res = fct::mult_pt(coeffs.at(step.key1), fpowers.at(step.key2), eb);
                    }
                    else if(step.key1[0] == 'P' && step.key2[0] == 'T')
                    {
                        res = fct::mult(fpowers.at(step.key1), fterms.at(step.key2), eb);
                    }
                    else if(step.key1[0] == 'T' && step.key2[0] == 'T')
                    {
                        res = fct::mult(fterms.at(step.key1), fterms.at(step.key2), eb);
                    }
                    else if(step.key1[0] == 'T' && step.key2[0] == 'P')
                    {
                        res = fct::mult(fterms.at(step.key1), fpowers.at(step.key2), eb);
                    }
                    else if(step.key1[0] == 'P' && step.key2[0] == 'C')
                    {
                        res = fct::mult_pt(coeffs.at(step.key2), fpowers.at(step.key1), eb, true);
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

    Ciphertext evaluate_cl(Ciphertext& x, ErrBound& eb)
    {
        Ciphertext fctx = x;
        
        // x^2
        Ciphertext x2 = square(fctx, eb);
        
        // -2x+3
        Ciphertext x3 = mult_pt(-2, x, eb, false);
        x3 = add_pt(x2, 3);

        // final
        Ciphertext res = mult(x2, x3, eb);
        return res;
    }
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

    if (e == 0.0) {
        return std::to_string(b);
    }

    const double log2e = std::log2(e);

    return std::format("{} {}2^({:.6f})", b, sign, log2e);
}