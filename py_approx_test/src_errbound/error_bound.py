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
    def cal_bound(self, m: float, coeff: list, err_mode: int) -> float:
        max_deg = len(coeff) - 1
        a = abs(coeff[max_deg])
        if err_mode == 1:
            if max_deg == 2:
                err = a * B_2(m, self.Bc, self.Bs) + self.Bs
            elif max_deg == 4:
                b = abs(coeff[max_deg - 2])
                err = (
                    (2 * a * pow(m, 3) * self.Bc)
                    + (2 * b * m * self.Bc)
                    + (a * pow(m, 2) * self.Bs)
                    + (2 * a * pow(m, 3) * self.Bc)
                    + ((a + 1) * pow(m, 2) * self.Bs)
                    + self.Bs
                )
            elif max_deg == 6:
                raise NotImplementedError("cal_bound case for degree 6 is not implemented.")
            else:
                raise ValueError("Unsupported maximum degree.")

        elif err_mode == 2:
            if max_deg == 2:
                err = self.Bs * (a + 1)
            elif max_deg == 4:
                B2 = self.Bs * (a + 1)
                b = abs(coeff[max_deg - 2])
                err = B2 * pow(m, 2) + self.Bs * (a * pow(m, 2) + b + 1)
                # err = (
                #     (a + 1) * self.Bs * pow(m, 2)
                #     # + (1 / self.scale) * (a + 1) * pow(self.Bs, 2)
                #     + self.Bs
                # )
            elif max_deg == 6:
                B2 = self.Bs * (a + 1)
                b = abs(coeff[max_deg - 2])
                B4 = B2 * pow(m, 2) + self.Bs * (a * pow(m, 2) + b + 1)
                c = abs(coeff[max_deg - 4])
                err = B4 * pow(m, 2) + (a * pow(m, 4) + b * pow(m, 2) + c + 1)
                # err = (
                #     (a + 1) * self.Bs * pow(m, 4)
                #     + pow(m, 2) * self.Bs
                #     + self.Bs
                # )
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
