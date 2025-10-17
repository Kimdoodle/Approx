from src_approx.approx_helut import helut
from src_approx.approx_remez import multi_remez
from src_approx.approx_sgn import sgn
from src_remez.print_plot import remove_images
from src_errbound.error_bound import EB
from parse import build_desmos_script


criteria = "ctxt" # depth, ctxt

###
# helut(p_num, e_num, interval_Scale)
# print("#" * 20)


# print("#" * 20)


# coeff = remez(p_num, e_num, interval_nonScale, True)
# print("#" * 20)

# remove_images()
max_n = 6
N = pow(2, 17)
sigma = 3.1
h = 192
s_num = 50
eb = EB(sigma, N, h, s_num)

for p_num in [5, 6, 7, 8, 9, 10]:
    for e_num in [30]:
        try:
            # sgn(p_num, e_num, interval_Scale, criteria)
            multi_remez(p_num, e_num, max_n, eb, "normal")
            print(build_desmos_script(f"doc/coeff_{p_num}_{e_num}.txt"))
        except Exception as e:
            # print("main")
            print(e)
            continue

# multi_remez(p_num, e_num, max_n, False, "normal")
# print(build_desmos_script(f"coeff_{p_num}_{e_num}.txt"))