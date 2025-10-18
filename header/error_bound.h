#pragma once
#include <cmath>
#include <vector>

class EB {
public:
    EB(double sigma, int N, int h)
        : Bc(B_clean(sigma, N, h)),
          Bs(B_scale(N, h)),
          B2(0.0), B4(0.0), B6(0.0) {}

    void cal(double m) {
        B2 = B_2(m, Bc, Bs);
        B4 = B_4(m, Bc, Bs);
        B6 = B_6(m, Bc, Bs);
    }

    // Error bound of coeff.
    double cal_bound(double m, const std::vector<double>& coeff) {
        cal(m);
        int max_deg = static_cast<int>(coeff.size()) - 1;
        double a = coeff[max_deg];
        double b = coeff[max_deg-2];
        switch(max_deg)
        {
            case 2:
                return 2 * a * m * Bc + (a+1) * Bs;
                break;
            case 4:
                return (2 * a * pow(m, 3) * Bc) + (2 * b * m * Bc) + (a * pow(m, 2) * Bs) + (2 * a * pow(m, 3) * Bc) + ((a+1) * pow(m, 2) * Bs) + Bs;
                break;
            case 6:
                return 0.0;
                break;
        }
    }
    
    double Bc, Bs, B2, B4, B6;

private:
    // Bclean
    static double B_clean(double sigma, int N, int h) {
        return (8 * sqrt(2.0) * sigma * N)
             + (6 * sigma * sqrt(static_cast<double>(N)))
             + (16 * sigma * sqrt(static_cast<double>(h) * N));
    }

    // Bscale
    static double B_scale(int N, int h) {
        return sqrt(N / 3.0) * (3.0 + 8.0 * sqrt(static_cast<double>(h)));
    }

    // B_2 = 2mB + 2B_s
    static double B_2(double m, double B_c, double B_s) {
        return (2.0 * m * B_c) + B_s;
        // return 2 * B_s;
    }

    // B_4 = 4m^3B + 4m^2B_s + 2B_s
    static double B_4(double m, double B_c, double B_s) {
        return (4.0 * pow(m, 2) * B_c) + (2.0 * m * B_s) + B_s;
        // return (4 * pow(m, 2) * B_s) + (2 * B_s);
    }

    // B_6 = 6m^5B + 6m^4B_s + 2m^2B_s + 2B_s
    static double B_6(double m, double B_c, double B_s) {
        return (4.0 * pow(m, 3) * B_c)
             + (2.0 * pow(m, 2) * B_c)
             + (2.0 * pow(m, 2) * B_s)
             + (2.0 * m * B_s)
             + B_s;
        // return (6 * pow(m, 4) * B_s) + (2 * pow(m, 2) * B_s) + (2 * B_s);
    }
};
