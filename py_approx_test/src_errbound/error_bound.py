from math import sqrt, pow, log2

class EB:
    def __init__(self, sigma: float, N: int, h: int, s: int):
        self.Bc = B_clean(sigma, N, h)
        self.Bs = B_scale(N, h)
        self.scale = pow(2, s)
        # print(log2(self.Bc/self.scale))
        # print(log2(self.Bs/self.scale))
        self.B2 = 0.0
        self.B4 = 0.0
        self.B6 = 0.0

    # def cal(self, m: float):
    #     print((2 * m * self.Bc + self.Bs) / self.scale)
    #     print(log2((2 * m * self.Bc + self.Bs) / self.scale))
    #     self.B2 = B_2(m, self.Bc, self.Bs)
    #     self.B4 = B_4(m, self.Bc, self.Bs)
    #     self.B6 = B_6(m, self.Bc, self.Bs)

    # Error bound of coeff.
    def cal_bound(self, m: float, coeff: list, coeff_type: str) -> float:
        max_deg = len(coeff) - 1
        a = abs(coeff[max_deg])
        
        match coeff_type:
            case "odd":
                if max_deg == 1:
                    err = self.Bs
                elif max_deg == 3:
                    B2 = self.Bs * (a + 1)
                    err = m * B2 + self.Bs
                elif max_deg == 5:
                    B2 = self.Bs * (a + 1)
                    B3 = m * B2 + self.Bs
                    err = m * B3 + self.Bs
            case "even":
                if max_deg == 2:
                    err = self.Bs * (a + 1)
                elif max_deg == 4:
                    B2 = self.Bs * (a + 1)
                    b = abs(coeff[max_deg - 2])
                    err = B2 * pow(m, 2) + self.Bs * (a * pow(m, 2) + b + 1)
                elif max_deg == 6:
                    B2 = self.Bs * (a + 1)
                    b = abs(coeff[max_deg - 2])
                    B4 = B2 * pow(m, 2) + self.Bs * (a * pow(m, 2) + b + 1)
                    c = abs(coeff[max_deg - 4])
                    err = B4 * pow(m, 2) + (a * pow(m, 4) + b * pow(m, 2) + c + 1)
            case "cl":
                b = abs(coeff[max_deg-1])
                return (a * m * self.Bs + self.Bs * self.Bs + b * self.Bs) / self.scale

            case "etc":
                match max_deg:
                    case 2:
                        err = (m + 1) * self.Bs
                    case 3:
                        a = abs(coeff[max_deg])
                        b = abs(coeff[max_deg - 1])
                        pt1, pt2 = pow(m, 2), (a*m + b)
                        er1, er2, er3 = self.Bs, self.Bs, self.Bs
                        err = (pt1 * er2) + (pt2 * er1) + er3 + self.Bs
                    case 4:
                        a = abs(coeff[max_deg])
                        b = abs(coeff[max_deg - 1])
                        c = abs(coeff[max_deg - 2])
                        pt1, pt2 = pow(m, 2), (a*pow(m, 2) + b*m + c)
                        er1, er2, er3 = self.Bs, (a+2) * self.Bs, self.Bs
                        err = (pt1 * er2) + (pt2 * er1) + er3 + self.Bs
                    case 5:
                        if a == 0:
                            raise ZeroDivisionError(f"a == 0  (coeff={coeff})")
                        a = abs(coeff[max_deg])
                        b = abs(coeff[max_deg - 1])
                        c = abs(coeff[max_deg - 2])
                        d = abs(coeff[max_deg - 3])
                        pt1, pt2 = a * pow(m, 3), (pow(m, 2) + (b/a * m) + c/a)
                        er1, er2 = (pow(m, 2) + (a * m) + 1) * self.Bs, 2 * self.Bs
                        er3 = (d + 2) * self.Bs
                        err = (pt1 * er2) + (pt2 * er1) + er3 + self.Bs
                    case 6:
                        a = abs(coeff[max_deg])
                        b = abs(coeff[max_deg - 1])
                        c = abs(coeff[max_deg - 2])
                        d = abs(coeff[max_deg - 3])
                        e = abs(coeff[max_deg - 4])
                        pt1, pt2 = (pow(m, 2) * (a*m + b) + (c*m + d)), pow(m, 3)
                        er1, er2, er3 = (pow(m, 2) + (a*m) + b + 2) * self.Bs, (m + 1) * self.Bs, (e + 2) * self.Bs
                        err = (pt1 * er2) + (pt2 * er1) + er3 + self.Bs
        return err
    
    def cal_bound_cleanse(self, m: float) -> float:
        return (2 * m + 4) * self.Bs

# Bclean
def B_clean(sigma: float, N: int, h: int) -> float:
    return (8 * sqrt(2.0) * sigma * N) + (6 * sigma * sqrt(N)) + (16 * sigma * sqrt(h * N))

# Bscale
def B_scale(N: int, h: int) -> float:
    return sqrt(N / 3.0) * (3.0 + 8.0 * sqrt(h))

# B_2 = 2mB + B_s
def B_2(m: float, B_c: float, B_s: float) -> float:
    return (2.0 * m * B_c) + B_s

# B_4 = 4m^2B + 2mB_s + B_s
def B_4(m: float, B_c: float, B_s: float) -> float:
    return (4.0 * pow(m, 2) * B_c) + (2.0 * m * B_s) + B_s

# B_6 = 4m^3B + 2m^2B + 2m^2B_s + 2mB_s + B_s
def B_6(m: float, B_c: float, B_s: float) -> float:
    return (4.0 * pow(m, 3) * B_c) + (2.0 * pow(m, 2) * B_c) + (2.0 * pow(m, 2) * B_s) + (2.0 * m * B_s) + B_s


