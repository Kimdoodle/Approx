#include "header/main.h"
#include <math.h>
#include "header/error_bound.h"

int main()
{
    int p_num = 2;
    int e_num = 30;
    int s_num = 50;
    int N_num = 17;
    size_t hwt = 192;
    double sigma = 3.1;
    double scale = pow(2, s_num);
    int N = pow(2, N_num);
    EB eb = EB(sigma, N, hwt);

    for(p_num=6; p_num<=7; p_num++)
        // test_base(p_num, s_num, e_num, N_num, hwt, true, false, true, true);
        test_time(p_num, s_num, e_num, N_num, false, true);
    // test_cleanse_s(p_num, e_num, s_num);

    // vector<double> coeff = {0.6804462163512842, 0.0, -0.4010128038041, 0.0, 0.040100790861333446};
    // test_err(eb, coeff, p_num, e_num, s_num, N_num, hwt);

}