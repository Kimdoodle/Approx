# from src_approx.approx_helut import helut
# from src_approx.approx_remez import multi_remez
# from src_approx.approx_sgn import sgn
from src_remez.algorithm import remez_algorithm
# from src_remez.print_plot import remove_images
from src_approx.approx_remez_MT import multi_remez
from src_errbound.error_bound import EB
from parse import build_desmos_script
import numpy as np
import os

criteria = "ctxt"
max_n = 6
N = pow(2, 17)
sigma = 3.1
h = 192
s_num = 50
eb = EB(sigma, N, h, s_num)


# killall -u doodle python
os.makedirs("doc", exist_ok=True)
for p_num in [2]:
    for e_num in [30]:
        try:
            multi_remez(p_num, e_num, max_n, eb, "normal", True)
            print(build_desmos_script(f"doc/coeff_{p_num}_{e_num}.txt"))
        except Exception as e:
            print("MAIN")
            print(e)
            continue
        # try:
        #     print(build_desmos_script(f"doc/coeff_{p_num}_{e_num}.txt"))
        # except:
        #     pass