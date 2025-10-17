import numpy as np
from src_remez.algorithm import remez_algorithm, cleanse
from src_remez.print_plot import *
from src_remez.print import print_intervals, debug_print, coeff2txt, interval2txt
from src_errbound.error_bound import EB
from cal_depth import cal_coeff, CalData
from math import log2
import copy

class remezData:
    def __init__(self, interval):
        self.coeff_log: list = []
        self.interval_log: list = [interval]
        self.CalData_log: list[CalData] = []
        self.pre_log: list[int] = []
        self.total_CalData = CalData()
        self.iter_cleanse = 0

    def update(self, coeff, pre, interval):
        self.add(coeff)
        self.pre_log.append(pre)
        self.interval_log.append(interval)

    def add(self, coeff: list):
        self.coeff_log.append(coeff)
        cal = cal_coeff(coeff)
        self.CalData_log.append(cal)
        self.total_CalData.add(cal, "add")

    def compare(self, challenger, eb: EB) -> int:
        if type(challenger) == list:
            challenger_cal = cal_coeff(challenger)
            res = self.total_CalData.compare(challenger_cal, "both")
            return res
    
        elif type(challenger) == remezData:
            res = self.total_CalData.compare(challenger.total_CalData, "both")
            if res == 0:
                if len(self.coeff_log) < len(challenger.coeff_log):
                    return 1
                elif len(self.coeff_log) > len(challenger.coeff_log):
                    return -1
                else:
                    if self.pre_log[-1] < challenger.pre_log[-1]:
                        return 1
                    elif self.pre_log[-1] > challenger.pre_log[-1]:
                        return -1
            return res

    """
        새로운 비교방안????
        1. 동일하게 각각 계산복잡도 비교
        2. 동일한 depth ->  
    """
    # def compare(self, challenger: 'remezData') -> int:
    #     pass

    def print_params(self):
        print(f"precision: ", end='')
        for pre in self.pre_log:
            print(f"{pre}", end='')
            if pre != self.pre_log[-1]:
                print(" -> ", end='')
            else:
                print()
        print(f"coeff: ", end='')
        for i, coeff in enumerate(self.coeff_log):
            print(f"{len(coeff)-1}", end='')
            if i+1 < len(self.coeff_log):
                print(" -> ", end='')
            else:
                print()
        print(f"total complexity:")
        self.total_CalData.print_params()
        print("#"*20)
    
def remez_recursion(routes: list, data: remezData, e_num: int, max_n: int, eb: EB, printmode="debug"):
    n = 3
    intervals = data.interval_log[-1]

    def evalF(x):
        return 1 if x >= 0.5 else 0

    while n <= max_n:
        if intervals[0][1] > intervals[1][0]:
            return
        coeff, max_err, next_intervals = remez_algorithm(n, intervals, evalF, "odd", eb, 2, "normal")

        pre = -int(log2(max_err))
        if coeff == [-1] or pre < 0:
            break
        
        debug_print(f"n={n}, pre {pre}, step {len(data.coeff_log) + 1}, intervals: ", printmode, "")
        if printmode == "debug":
            print_intervals(next_intervals, 10)

        # 새로운 remezData 객체 생성 후 정보 갱신
        newData = copy.deepcopy(data)
        newData.update(coeff, pre, next_intervals)

        if pre >= e_num: # 성공
            routes.append(newData)
            return
        
        # Cleanse 고려 - 최대 2번
        newData_cl = newData
        for iter in range(2):
            max_err_cl, next_intervals_cl = cleanse(eb, next_intervals)
            pre_cl = -int(log2(max_err_cl))
            if pre >= pre_cl:
                break
            else:
                max_err = max_err_cl
                pre = pre_cl
                next_intervals = next_intervals_cl

            debug_print(f"n={n}, pre {pre}, step {len(data.coeff_log) + 1} +  cleanse {iter+1}\n \t\tintervals: ", printmode, "")
            if printmode == "debug":
                print_intervals(next_intervals, 10)

            newData_cl = copy.deepcopy(newData_cl)
            newData_cl.update([0, 0, 3, -2], pre, next_intervals)
            newData_cl.iter_cleanse += 1
            if pre >= e_num:
                routes.append(newData_cl)
                return
                
        remez_recursion(routes, newData, e_num, max_n, eb, printmode) # 새 객체가 정확도 조건을 만족할 때까지 재귀

        n += 2
    return

def multi_remez(p_num: int, e_num: int,max_n: int, eb: EB, printmode: str):    
    p = pow(2, p_num)
    e = eb.Bc/eb.scale    
    intervals = [[0, e], [(1-e), (p-1+e)]]   
    
    routes = [] # 최종 경로 목록

    # 초기 evalF
    def evalF(x):
        return 1 if x < 0.5 else 0
        
    # 첫구간 remez 수행 후 재귀
    for n in range(2, max_n+1, 2):
        coeff, max_err, next_intervals = remez_algorithm(n, intervals, evalF, "even", eb, 2, "normal")
        pre = -int(log2(max_err))
        if max_err == 99:
            debug_print(f"Init approximation n={n} failed.", printmode)
            continue
        debug_print(f"Init approximation\nn={n}, step 0, intervals: ", printmode, "")
        if printmode == "debug":
            print_intervals(next_intervals, 10)

        # remezData 객체 생성 후 정보 갱신
        data = remezData(intervals)
        data.update(coeff, pre, next_intervals)

        # 재귀 - 홀수차로 진행
        remez_recursion(routes, data, e_num, max_n, eb, printmode)
    
    if len(routes) > 0 and routes[0].pre_log[-1] >= e_num: # 성공
        print(f"Total {len(routes)} routes.")

        # 최종 비교
        bob = [routes[0]]

        for i in range(1, len(routes)):
            route = routes[i]
            res = bob[0].compare(route, eb)
            if res == -1:
                bob = [route]
            elif res == 0:
                bob.append(route)

        # 텍스트파일 저장
        with open(f"doc/coeff_{p_num}_{e_num}.txt", "w", encoding="utf-8") as f:
            for coeff in bob[-1].coeff_log:
                f.write(coeff2txt(coeff))
                f.write('\n')
        
        print(f"(p, e): {p_num}, {e_num}")
        # for b in bob:
            # b.print_params()
        bob[-1].print_params()

    else: # 실패
        print("Approximation Failed.")
        # routes[0].print_params()