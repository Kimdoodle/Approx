from src_approx.approx_remez_MT import multi_remez
from src_approx.util_MT import set_max_workers
from src_errbound.error_bound import EB
from parse import build_desmos_script
import os

criteria = "ctxt"
max_n = 6
N = pow(2, 17)
sigma = 3.1
h = 192
s_num = 50
eb = EB(sigma, N, h, s_num)
tail = "oe"
printmode = "normal"
multiThread = True
useDP = True

# killall -u doodle python

set_max_workers(32)
os.makedirs("doc", exist_ok=True)
for p_num in range(9, 10):
    for e_num in [30]:
        for _ in range(5):
            try:
                multi_remez(p_num, e_num, max_n, eb, tail, printmode, multiThread, useDP)
                print(build_desmos_script(f"doc/coeff_{tail}_{p_num}_{e_num}.txt"))
                break
            except Exception as e:
                print("MAIN")
                print(e)
                continue
        # try:
        #     print(build_desmos_script(f"doc/coeff_{p_num}_{e_num}.txt"))
        # except:
        #     pass