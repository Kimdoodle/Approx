from itertools import product
import math
from cal_depth import *
from src_remez.basic import evalP

# Sgn functions
f1 = [0, 3/2, 0, -1/2]
f2 = [0, 15/8, 0, -10/8, 0, 3/8]
f3 = [0, 35/16, 0, -35/16, 0, 21/16, 0, -5/16]
f4 = [0, 315/128, 0, -420/128, 0, 378/128, 0, -180/128, 0, 35/128]
coeffs = {
    1: f1, 
    2: f2, 
    3: f3,
    4: f4
}

def iter_fk(f_num: int, iter: int, x: float) -> float:
    y = x
    for _ in range(iter):
        y = evalP(coeffs[f_num], y)
    return y

def iter_sgn(f_num: int, e: float, iter: int, x: float) -> float:
    data1 = iter_fk(f_num, iter, x+2*e)
    data2 = iter_fk(f_num, iter, x-2*e)
    return 0.5*(data1 - data2)

def find_best(data_list: list[CalData], criteria: str):
    if not data_list:
        return []

    best_list = [data_list[0]]

    for candidate in data_list[1:]:
        comparison = candidate.compare(best_list[0], criteria)
        if comparison == 1:
            best_list = [candidate]
        elif comparison == 0:
            best_list.append(candidate)

    return best_list

def sgn(p_num: int, e_num: int, intervals: list, criteria: str):
    p = pow(2, p_num)
    e = pow(2, -e_num)
    f_nums = (1, 2, 3, 4)

    # 정밀도 조건을 만족하는 최소 반복횟수 검색
    results: list[CalData] = []
    for f_num in f_nums:
        for iter in range(1, 100):
            # 구간 [0, e/p]
            v0 = iter_sgn(f_num, e/p, iter, 0)
            v1 = iter_sgn(f_num, e/p, iter, intervals[1][1])
            cond_left  = (v0 >= 1-e) and (v1 >= 1-e)

            # 구간 [(1-e)/p, 1]
            v2 = iter_sgn(f_num, e/p, iter, intervals[2][0])
            v3 = iter_sgn(f_num, e/p, iter, intervals[2][1])
            if v2 >= 0.5 or v3 >= 0.5:
                pass
                # print(f"Error! v2={v2}, v3={v3}")
            cond_right = (v2 <= e) and (v3 <= e)    

            if cond_left and cond_right:
                result = cal_iter(coeffs[f_num], iter)
                result.title = f"f_{f_num}(x)"
                result.iter = iter
                results.append(result)
                break
    
    # f_1 ~ f_4 중 최적의 결과 출력
    best_results = find_best(results, criteria)
    
    print("Sgn+Sgn Result")
    for data in best_results:
        # 텍스트파일 저장
        with open(f"sgn.txt", "a", encoding="utf-8") as f:
            f.write(f"p={p_num}\t\te={e_num}\t\t{data.title}\t\t{data.iter}")
            f.write('\n')
        # print(f"f_{i+1}(x): iter {data[0]} times.")
        data.print_params()
        print('-'*10)




if __name__ == '__main__':
    p = 4
    e = 3
    sgn(p, e)
