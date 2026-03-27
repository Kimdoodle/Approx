// decomp.cpp
#include "../include/polyEval_class.h"

// Class Complexity
Complexity::Complexity() = default;

void Complexity::insert_value(int depth, int cmult, int pmult, int add)
{
    this->depth = depth;
    this->cmult = cmult;
    this->pmult = pmult;
    this->add = add;
}

// Class Poly
Poly::Poly(std::vector<double> coeff)
{
    this->coeff = trim_poly(coeff);
    this->deg = std::max((int)this->coeff.size() - 1, 0);
    this->complexity = Complexity();
}

std::pair<Poly, Poly> Poly::seperate_poly(int i, bool multA)
{
    std::vector<double> coeff_p, coeff_q;

    // slicing, assuming value i won't cause any errors.
    coeff_p.assign(this->coeff.begin() + i, this->coeff.end());
    coeff_q.assign(this->coeff.begin(), this->coeff.begin() + i);

    // trim
    coeff_p = trim_poly(coeff_p);
    coeff_q = trim_poly(coeff_q);

    // multA
    if (multA && !coeff_p.empty()) {
        double leading = coeff_p.back();
        for (double &c : coeff_p) {
            c /= leading;
        }
    }

    return {Poly(coeff_p), Poly(coeff_q)};
}

// Class XI
XI::XI(bool multA, int n)
{
    this->multA = multA;
    this->n = n;
    this->made_powers = {0, 1};
    int val = multA ? 1 : 0;
    this->depth = (int)ceil(log2(n + val));
    this->pmult = val;
}

void XI::add_routes(std::vector<std::pair<int, int>> route, std::set<int> made_powers)
{
    this->route = route;
    this->add_count = (int)route.size();    
    this->made_powers = made_powers;
}

// Class Deocmp
Decomp::Decomp(std::vector<double>& c, Complexity comp, XI xi) 
{
    coeff = c;
    // this->comp = comp;
    this->xi = xi;
    this->made_powers = this->xi.made_powers;
}

void Decomp::update(XI xi, std::shared_ptr<Decomp> dcmp_p, std::shared_ptr<Decomp> dcmp_q)
{
    this->xi = xi;
    this->dcmp_p = dcmp_p;
    this->dcmp_q = dcmp_q;
    this->made_powers.insert(xi.made_powers.begin(), xi.made_powers.end());
}

std::string Decomp::restore_dcmp() {
    if (this->coeff.empty()) 
        return "";

    // i=0
    if (this->xi.n == 0) {
        return vec_to_str(this->coeff); 
    }
    else
    {
        std::string str_i = "";
        std::string str_p = "";
        std::string str_q = "";
        double leading_coeff = 1.0;
        std::vector<double> coeff_p = this->dcmp_p->coeff;
        
        // polynomial seperation.
        int coeff_p_max_index = (int)coeff_p.size() - 1;
        if(this->xi.multA && coeff_p.size() > 0)
        {
            leading_coeff = coeff_p[coeff_p_max_index];
            for(int i=0; i<=coeff_p_max_index; i++)
                coeff_p[i] /= leading_coeff;
        }
        
        // restore ax^i
        std::string coeff_i = !this->xi.multA ? "" : std::to_string(coeff.back());
        if(this->xi.n == 1)
            str_i += ("(" + coeff_i + "x)");
        else
            str_i += ("(" + coeff_i + "x^" + std::to_string(this->xi.n) + ")");
        
        // restore p(x)
        if(!coeff_p.empty())
            str_p = "[ " + this->dcmp_p->restore_dcmp() + " ]";

        // restore q(x)
        if(this->dcmp_q)
            str_q = " + (" + this->dcmp_q->restore_dcmp() + ")";

        return str_i + str_p + str_q;
    }
}

void Decomp::merge_mp()
{
    if (dcmp_p) {
        dcmp_p->merge_mp();
        made_powers.insert(dcmp_p->made_powers.begin(), dcmp_p->made_powers.end());
    }

    if (dcmp_q) {
        dcmp_q->merge_mp();
        made_powers.insert(dcmp_q->made_powers.begin(), dcmp_q->made_powers.end());
    }
}

std::vector<double> trim_poly(const std::vector<double>& coeff) {
    std::vector<double> res = coeff;
    while (!res.empty() && res.back() == 0.0) {
        res.pop_back();
    }
    return res;
}

static bool is_integer_like(double x) {
    return std::isfinite(x) && x == std::trunc(x);
}

static std::string format_coeff(double x) {
    std::ostringstream oss;
    if (is_integer_like(x)) {
        // 정수(소수점 제거)
        oss << static_cast<long long>(std::llround(x));
    } else {
        // 소수
        oss << x;
    }
    return oss.str();
}

// std::string vec_to_str(const std::vector<double>& poly) {
//     std::string res;

//     if (poly.empty()) return res;

//     for (int i = static_cast<int>(poly.size()) - 1; i >= 0; --i) {
//         double coeff = poly[static_cast<size_t>(i)];
//         if (coeff == 0.0) continue;

//         std::string mark;
//         if (i == static_cast<int>(poly.size()) - 1 || coeff < 0.0) {
//             mark = "";
//         } else {
//             mark = "+";
//         }

//         std::string c = format_coeff(coeff);

//         if (i == 0) {
//             res += mark + c;
//         } else if (i == 1) {
//             res += mark + c + "x";
//         } else {
//             res += mark + c + "(x^" + std::to_string(i) + ")";
//         }
//     }
//     return res;
// }

std::string vec_to_str(const std::vector<double>& poly) {
    std::string res;
    if (poly.empty()) return res;

    for (int i = static_cast<int>(poly.size()) - 1; i >= 0; --i) {
        double coeff = poly[static_cast<size_t>(i)];
        if (coeff == 0.0) continue;

        // sign/mark
        std::string mark;
        if (res.empty()) {
            mark = (coeff < 0.0) ? "-" : "";
        } else {
            mark = (coeff < 0.0) ? "-" : "+";
        }

        // abs(coeff) string for display
        const double abs_coeff = std::abs(coeff);

        // 계수가 1인 경우(i>0) 계수 표시 생략 (단, -1은 '-'만 남도록)
        std::string c;
        if (i > 0 && is_integer_like(abs_coeff) && std::llround(abs_coeff) == 1) {
            c = "";  // 1x, 1(x^n) -> x, (x^n)
        } else {
            c = format_coeff(abs_coeff);
        }

        if (i == 0) {
            // 상수항은 1이라도 표시해야 함: 1, -1
            res += mark + format_coeff(abs_coeff);
        } else if (i == 1) {
            res += mark + c + "x";
        } else {
            res += mark + c + "(x^" + std::to_string(i) + ")";
        }
    }
    return res;
}

Complexity attach_poly(Complexity c1, Poly d2, Complexity c2, char attach_type)
{
    /*
    attach
    두 Complexity를 attach_type에 따라 결합
    attach_type = 'x' --> 곱셈. depth는 max(c1.d, c2.d) + 1, 나머지는 단순 덧셈.
    attach_type = '+' --> 덧셈. depth는 max(c1.d, c2.d) 이며 cmult, pmult는 단순 덧셈.
                                + add는 (덧셈 + 1)
    */
    Complexity res = Complexity();
    switch(attach_type)
    {
    case '+':
        if(d2.deg == 0 && d2.coeff[0] == 0.0)
            return c1;
        res.insert_value(std::max(c1.depth, c2.depth), c1.cmult + c2.cmult, c1.pmult + c2.pmult, c1.add + c2.add + 1);
        break;
    
    case 'x':
        res.insert_value(std::max(c1.depth, c2.depth) + 1, c1.cmult + c2.cmult, c1.pmult + c2.pmult, c1.add + c2.add);
        break;
    }
    return res;
}

Complexity mult_xi(XI d1, Poly d2, Complexity c2)
{
    /*
    mult - XI클래스와 Poly의 곱셈
    depth는 max(c1.d, c2.d) + 1 이며 나머지는 단순 덧셈
    */
    Complexity res = Complexity();
    res.insert_value(std::max(d1.depth, c2.depth) + 1, d1.add_count + c2.cmult, d1.pmult + c2.pmult, c2.add);
    return res;
}