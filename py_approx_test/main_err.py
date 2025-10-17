# from src_errbound.error_bound import EB

# N = pow(2, 17)
# h = 192
# sigma = 3.6
# x_bound = pow(2, 6)

# coeff = [0, 0, 1, 0, 1]

# eb = EB(sigma, N, h)

# print(eb.cal_bound(1, coeff))

Bc=4.309954427250231e-09
Bs=2.1136465360237085e-11

from math import log2

B0 = 4*Bc+Bs
print(B0)