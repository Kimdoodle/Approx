from src_remez.basic import evalP
from cal_depth import cal_iter, CalData
from math import log2

def sqmethod(p: int, u: float, e: float) -> tuple[int, int, tuple, CalData]:
    iter = 0
    b1 = 1 - 2 * (e**2) / (p**2)
    b2 = 1 - 2 * ((1-e)**2) / (p**2)
    b3 = 1 - 2 * ((p-1+e)**2) / (p**2)

    while True:
        if (b1 > 1-u) and (b2 < u) and (b3 < u):
            break
        
        b1 = b1 * b1
        b2 = b2 * b2
        b3 = b3 * b3
        iter += 1
        # print(f"{b1}\t{b2}\t{b3}")
        
        err = max(1-b1, b2, b3)
        u_num = -int(log2(err))


    # 연산량 계산
    cal_res = CalData()
    cal_res.add_count = 1
    cal_res.pt_mult_count = 1
    cal_res.ct_mult_count = iter + 1
    cal_res.depth = iter + 2
    
    return iter, u_num, (b1, b2, b3), cal_res

def cleanse(e: float, data: tuple) -> tuple:
    coeff = [0.0, 0.0, 3.0, -2.0]
    d1, d2, d3 = data
    if d1 >= 1-e and d2 <= e and d3 <= e:
        return 0, CalData()
    s = 0
    
    while True:
        d1 = evalP(coeff, d1)
        d2 = evalP(coeff, d2)
        d3 = evalP(coeff, d3)
        s += 1
        # print(f"{res1}\t{res2}\t{res3}")
        if d1 >= 1 - e and d2 <= e and d3 <= e:
            break

    # 연산량 계산
    cal_res = cal_iter(coeff, s)
    return s, cal_res

def helut(p_num: int, e_num: int):
    p = pow(2, p_num)
    e = pow(2, -e_num)

    best = (0, 0, 0) # r, u, s
    best_caldata = CalData()
    best_caldata.depth = 99
    best_caldata.add_count = 99
    best_caldata.pt_mult_count = 99
    best_caldata.ct_mult_count = 99

    # 모든 경우의 수에 대하여 계산
    for u_num in range(1, e_num):
        # print(f"u : 2^-{u_num}")
        u = pow(2, -u_num)
        r, u_num, data, cal = sqmethod(p, u, e)
        s, cal2 = cleanse(e, data)
        # print(f"r, s:\t({r}, {s})")

        cal.add(cal2)
        # cal.print_params()
        # print('-'*10)
        if s > 0:
            res = best_caldata.compare(cal, "both")

            if res < 0:
                best = (r, u_num, s)
                best_caldata = cal

    # 결과 출력
    print("HELUT result")
    print(f"r, u, s:\t({best[0]}, {best[1]}, {best[2]})")
    best_caldata.print_params()
    

if __name__ == '__main__':
    p_num = 5
    e_num = 30
    helut(p_num, e_num)


    