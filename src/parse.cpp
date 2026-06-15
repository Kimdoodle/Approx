// parse.cpp
#include "../include/parse.h"

static double stod_safe(const std::string& s) {
    // 필요시 안전 변환용 헬퍼
    return std::stod(s);
}

std::vector<std::vector<double>> parse_remez_coeff(int p_num, int e_num, std::string criteria, bool print)
{
    int line_count = 0;  // 전체 줄 개수

    // 1) 파일명: coeff_pnum_enum.txt 
    const std::string filename = std::format("doc/coeff_{}_{}_{}.txt", criteria, p_num, e_num);
    // const std::string filename = std::format("doc_ERE3/data/coeff_{}_{}_{}.txt", criteria, p_num, e_num);
    std::cout << std::format("Opening {}", filename) << std::endl;
    
    // 2) 계수 파일 파싱
    std::ifstream infile(filename);
    if (!infile)
        throw std::runtime_error("파일을 열 수 없습니다: " + filename);

    std::vector<std::vector<double>> result;
    std::string line;

    // 계수 추출 패턴: (계수)x^(지수)
    const std::regex term_re(
        R"(([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?))" // coeff (그룹1)
        R"(x\^)"                                              // x^
        R"((\d+))"                                            // exponent (그룹2)
    );

    // s 계산용 패턴: "+3x^2-2x^3" (공백 가변 허용, 앞 '+'는 선택)
    const std::regex s_re(R"((?:\+)?\s*3\s*x\^2\s*-\s*2\s*x\^3)");

    // 파일 순회
    while (std::getline(infile, line)) {
        ++line_count;  // 줄 개수 증가

        // (exp -> coeff) 맵 구성
        std::map<int, double> exp_to_coeff;
        int max_exp = -1;

        for (std::sregex_iterator it(line.begin(), line.end(), term_re);
             it != std::sregex_iterator(); ++it)
        {
            const std::string coeff_str = (*it)[1].str();
            const std::string exp_str   = (*it)[2].str();

            double c = stod_safe(coeff_str);
            int    e = std::stoi(exp_str);

            exp_to_coeff[e] = c;
            if (e > max_exp) max_exp = e;
        }

        // 항이 하나도 없으면 빈 벡터
        if (max_exp < 0) {
            result.emplace_back();
            continue;
        }

        // 누락 차수 0.0으로 채움
        std::vector<double> coeffs(static_cast<size_t>(max_exp + 1), 0.0);
        for (const auto& kv : exp_to_coeff) {
            coeffs[static_cast<size_t>(kv.first)] = kv.second;
        }

        result.push_back(std::move(coeffs));
    }
    return result;
}

int calcuate_remez_depth(std::vector<std::vector<double>> coeffs)
{
    int total = 0;
    for(auto coeff: coeffs)
    {
        total += ceil(log2(coeff.size()));
    }
    return total;
}

nlohmann::json load_decomp_cache(const std::string& filename)
{
    std::ifstream f(filename);
    
    // 파일 열기 실패 확인
    if (!f.is_open()) {
        std::cerr << "Failed to open " << filename << std::endl;
        return nullptr; // 실패 시 null 상태의 json 반환
    }

    nlohmann::json cache_root;
    try {
        cache_root = nlohmann::json::parse(f);
    } catch (nlohmann::json::parse_error& e) {
        std::cerr << "JSON Parse Error: " << e.what() << std::endl;
        return nullptr;
    }

    return cache_root;
}

std::string get_poly_type_key(const std::vector<double>& coeff) {
    std::string key;
    key.reserve(coeff.size());

    for (double v : coeff) {
        double ip;
        if(v == 0.0)
            key.push_back('0');
        else if(std::modf(v, &ip) == 0.0)
            key.push_back('I');
        else
            key.push_back('F');

        // key.push_back(v == 0.0 ? '0' : 'F');
    }
    return key;
}

std::shared_ptr<Decomp> reconstruct_decomp_from_cache(const json& plan, Poly poly, int depth) {
    // 1. Base Case: Check leaf node
    if (plan["is_leaf"].get<bool>())
        return std::make_shared<Decomp>(poly.coeff, Complexity());

    // 2. XI data 
    json plan_xi = plan["xi"];
    XI xi = XI(plan_xi.value("multA", false), plan_xi.value("n", 0));
    
    // made_powers 복원
    std::set<int> made_powers = {0, 1};
    std::vector<std::pair<int, int>> routes;
    for(auto& r : plan_xi.value("route", json::array())) {
        routes.push_back({r[0], r[1]});
    }    
    for(int i=0; i<routes.size(); i++)
    {
        if(i == routes.size()-1 && xi.multA)
        {
            break;
        }
        std::pair<int, int> op = routes[i];
        made_powers.insert(op.first + op.second);
    }
    xi.coeff = poly.coeff.back();
    xi.add_routes(routes, made_powers);
    
    // 3. Seperate polynomial
    std::pair<Poly, Poly> polys = poly.seperate_poly(xi.n, xi.multA);
    Poly poly_p = polys.first;
    Poly poly_q = polys.second;

    // 4. Reconstruct p(x), q(x)
    std::shared_ptr<Decomp> dcmp_p = reconstruct_decomp_from_cache(plan["p"], poly_p, depth+1);
    std::shared_ptr<Decomp> dcmp_q;
    if(!poly_q.coeff.empty())
        dcmp_q = reconstruct_decomp_from_cache(plan["q"], poly_q, depth+1);

    // 5. Complexity
    // Complexity comp_i = Complexity();
    // comp_i.insert_value(xi.depth, xi.add_count, xi.pmult, 0);
    // Complexity comp_p = dcmp_p->comp;
    // Complexity comp_pi = mult_xi(xi, poly_p, comp_p);
    // Complexity comp_piq = comp_pi;
    // if(dcmp_q)
    //     comp_piq = attach_poly(comp_pi, poly_q, dcmp_q->comp, '+');

    // 6. final result
    // std::shared_ptr<Decomp> res = std::make_shared<Decomp>(poly.coeff, comp_piq);
    std::shared_ptr<Decomp> res = std::make_shared<Decomp>(poly.coeff, Complexity());
    res->total_routes = xi.route;
    res->update(xi, dcmp_p, dcmp_q);

    // routes, made_powers 병합
    res->made_powers.insert(dcmp_p->made_powers.begin(), dcmp_p->made_powers.end());
    res->merge_route(dcmp_p->total_routes);
    if(dcmp_q)
    {
        res->made_powers.insert(dcmp_q->made_powers.begin(), dcmp_q->made_powers.end());
        res->merge_route(dcmp_q->total_routes);
    }
    return res;
}

void Decomp::merge_route(const std::vector<std::pair<int, int>>& a)
{
    std::set<std::pair<int, int>> res;
    res.insert(this->total_routes.begin(), this->total_routes.end());
    res.insert(a.begin(), a.end());

    this->total_routes = std::vector<std::pair<int, int>>(res.begin(), res.end());
}

// Class EvalStep
EvalStep::EvalStep() = default;

EvalStep::EvalStep(std::shared_ptr<Decomp> dcmp)
{
    this->add_value('P', 1, "x");
    //x^i들 구성
    for(std::pair<int, int> route: dcmp->total_routes)
    {
        if(dcmp->made_powers.contains(route.first + route.second))
        {
            std::string key1 = "P" + std::to_string(route.first);
            std::string key2 = "P" + std::to_string(route.second);
            this->add_step(key1, key2, 'x');
        }
    }
    serialize_dcmp(dcmp);
}

std::string EvalStep::serialize_dcmp(std::shared_ptr<Decomp> dcmp)
{    
    XI xi = dcmp->xi;
    std::string coeff_key, xi_key, axi_key;

    // leaf polynomial
    if(!dcmp->dcmp_p && !dcmp->dcmp_q)
    {
        std::string res_key = "";
        for(int i=dcmp->coeff.size()-1; i>=0; i--)
        {
            if(i == 0 && dcmp->coeff[i] != 0)
            {
                coeff_key = this->add_value('C', this->coeff_count++, std::format("{}", dcmp->coeff[0]));
                res_key = res_key == "" ? coeff_key : this->add_step(res_key, coeff_key, '+');
            }
            if(i != 0 && dcmp->coeff[i] != 0)
            {
                if(dcmp->coeff[i] != 1)
                {
                    coeff_key = this->add_value('C', this->coeff_count++, std::format("{}", dcmp->coeff[i]));
                    xi_key = std::format("P{}", i);
                    axi_key = this->add_step(coeff_key, xi_key, 'x');
                    res_key = res_key == "" ? axi_key : this->add_step(res_key, axi_key, '+');
                }
                else
                {
                    xi_key = std::format("P{}", i);
                    res_key = res_key == "" ? xi_key : this->add_step(res_key, xi_key, '+');
                }
            }
        }
        return res_key;
    }

    // multA여부에 따라 연산 순서 등록
    if(xi.multA)
    {
        char key_header = 'C';
        int key_number = this->coeff_count++;
        coeff_key = this->add_value(key_header, key_number, std::format("{}", xi.coeff));
        if(dcmp->made_powers.contains(xi.n))
        {
            //ax^n인 경우
            xi_key = std::format("P{}", xi.n);
            axi_key = this->add_step(coeff_key, xi_key, 'x');
        }
        else
        {
            if (xi.route.empty()) // a * x^i
            {
                xi_key = std::format("P{}", xi.n);
                axi_key = this->add_step(coeff_key, xi_key, 'x');
            }
            else // ax^i * x^j
            {
                xi_key = std::format("P{}", xi.route.back().first);
                axi_key = this->add_step(coeff_key, xi_key, 'x');
                std::string xj_key = std::format("P{}", xi.route.back().second);
                axi_key = this->add_step(axi_key, xj_key, 'x');
            }
        }
    } 
    else
        axi_key = std::format("P{}", xi.n);

    std::string axipx_key = axi_key;
    if(dcmp->dcmp_p)
    {
        std::string px_key = serialize_dcmp(dcmp->dcmp_p);
        axipx_key = this->add_step(axi_key, px_key, 'x');
    }

    std::string res_key = axipx_key;
    if(dcmp->dcmp_q)
    {
        std::string qx_key = serialize_dcmp(dcmp->dcmp_q);
        res_key = this->add_step(axipx_key, qx_key, '+');
    }
    return res_key;
}

std::string EvalStep::get_value(std::string key)
{
    // key "P{i}" --> returns data from powers.
    // key "C{i}" --> returns data from coeffs.
    // key "T{i}" --> returns data from terms.
    switch(key[0])
    {
        case 'P':
            return powers.at(key);
        case 'C':
            return coeffs.at(key);
        case 'T':
            return terms.at(key);
        default:
            throw std::invalid_argument("EvalStep::get_value: invalid key prefix: " + key);
    }
}

std::string EvalStep::add_value(char key_header, int key_number, std::string data)
{
    std::string key = key_header + std::to_string(key_number);
    switch(key_header)
    {
        case 'C':
            this->coeffs.insert_or_assign(key, data);
            this->eval_step.push_back(step{'o', key, data, ""});
            break;
        case 'P':
            this->powers.insert_or_assign(key, data);
            if(key_number == 1)
                this->eval_step.push_back(step{'o', key, data, ""});
            break;
        case 'T':
            this->terms.insert_or_assign(key, data);
            break;
    }
    return key;
}

std::string EvalStep::add_step(std::string key1, std::string key2, char op)
{
    char key_header;
    int term_number;
    std::string data, save_key;

    // 덧셈의 경우
    if(op == '+')
    {
        key_header = 'T';
        term_number = this->term_count++;
        data = std::format("({} + {})", this->get_value(key1), this->get_value(key2));
    }
    // x^i * x^j의 경우
    else if(key1[0] == 'P' && key2[0] == 'P')
    {
        key_header = 'P';
        term_number = (key1[1]-'0') + (key2[1]-'0');
        data = std::format("x^{}", term_number);
    }
    // a * x^i의 경우
    else if(key1[0] == 'C' && key2[0] == 'P')
    {
        key_header = 'T';
        term_number = this->term_count++;
        data = this->get_value(key1) + this->get_value(key2);
    }
    // x^i * f(x)의 경우
    else if(key1[0] == 'P' && key2[0] == 'T')
    {
        key_header = 'T';
        term_number = this->term_count++;
        data = this->get_value(key1) + this->get_value(key2);
    }
    // x^i * (f(x)=a)의 경우
    else if(key1[0] == 'P' && key2[0] == 'C')
    {
        key_header = 'T';
        term_number = this->term_count++;
        data = this->get_value(key1) + this->get_value(key2);
    }
    // ax^i * f(x)의 경우
    else if(key1[0] == 'T' && key2[0] == 'T')
    {
        key_header = 'T';
        term_number = this->term_count++;
        data = this->get_value(key1) + this->get_value(key2);
    }
    else if(key1[0] == 'T' && key2[0] == 'P')
    {
        key_header = 'T';
        term_number = this->term_count++;
        data = this->get_value(key1) + this->get_value(key2);
    }
    else
    {
        throw std::invalid_argument(std::format("EvalStep::add_step: invalid key data {} {} {}", key1, key2, op));
    }
    
    save_key = this->add_value(key_header, term_number, data);
    this->eval_step.push_back(step{op, key1, key2, save_key});
    return save_key;
}

void EvalStep::print_step()
{
    int count = 1;
    for(step s: this->eval_step)
    {
        this->print_line(s, count);
        ++count;
    }
}

void EvalStep::print_line(step s, int& count)
{
    std::string op_string;
    switch(s.op)
    {
        case '+':
            op_string = "ADD";
            // std::cout << std::format("{}\t{}\t{}\t{}\t{}\t{}\n", count++, op_string, s.key1, s.key2, s.save_key, this->get_value(s.save_key));
            break;
        case 'x':
            op_string = "MUL";
            // std::cout << std::format("{}\t{}\t{}\t{}\t{}\t{}\n", count++, op_string, s.key1, s.key2, s.save_key, this->get_value(s.save_key));
            break;
        case 'o':
            op_string = "RES";
            break;
    }
    std::cout << std::format("{}\t{}\t{}\t{}\t{}\n", count, op_string, s.key1, s.key2, s.save_key);
}

// std::variant<seal::Plaintext, seal::Ciphertext> get_value2(std::string key, std::map<int, seal::Ciphertext>& powers, std::map<int, seal::Plaintext>& coeffs, std::map<int, seal::Ciphertext>& terms)
// {

// }

// std::variant<seal::Plaintext, seal::Ciphertext> operate(std::string key1, std::string key2, char op, std::map<int, seal::Ciphertext>& powers, std::map<int, seal::Plaintext>& coeffs, std::map<int, seal::Ciphertext>& terms)
// {
//     char key1_header = key1[0];
//     char key2_header = key2[0];
//     int key1_number = key1[1];
//     int key2_number = key2[1];
// // 
//     if(key1_header == 'P')
//     {

//     }
    
// }