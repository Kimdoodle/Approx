import numpy as np
from src_remez.algorithm import remez_algorithm
from src_remez.print_plot import *
from src_remez.print import print_intervals, debug_print, coeff2txt, interval2txt
from src_remez.basic import calculate_next_remez
from src_errbound.error_bound import EB
from cal_depth import cal_coeff, CalData, check_coeff_type
from math import log2
from multiprocessing import Pool, cpu_count
import copy

class remezData:
    def __init__(self, interval):
        self.coeff_log: list = []
        self.interval_log: list = [interval]
        self.CalData_log: list[CalData] = []
        self.pre_log: list[tuple[float, float]] = []
        self.total_CalData = CalData()
        self.iter_cleanse = 0

    def update(self, coeff, pre: tuple, interval):
        self.add(coeff)
        self.pre_log.append(pre)
        self.interval_log.append(interval)

    def add(self, coeff: list):
        self.coeff_log.append(coeff)
        cal = cal_coeff(coeff)
        self.CalData_log.append(cal)
        self.total_CalData.add(cal, "add")

    def compare(self, challenger: 'remezData', mode="end") -> bool | int:
        self_depth = self.total_CalData.depth
        self_cmult = self.total_CalData.cmult
        self_pmult = self.total_CalData.pmult
        self_add = self.total_CalData.cadd
        chal_depth = challenger.total_CalData.depth
        chal_cmult = challenger.total_CalData.cmult
        chal_pmult = challenger.total_CalData.pmult
        chal_add = challenger.total_CalData.cadd

        self_comp = (self_depth, self_cmult, self_pmult, self_add)
        chal_comp = (chal_depth, chal_cmult, chal_pmult, chal_add)

        if self_comp > chal_comp:
            return True
        elif self_comp == chal_comp:
            return 99
        else:
            return False

    def print_params(self):
        print(f"precision: ", end='')
        for pre in self.pre_log:
            print(f"[{int(pre[0])}, {int(pre[1])}]", end='')
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
    
    def print_coeff_type(self):
        for coeff in self.coeff_log:
            ct = check_coeff_type(coeff)
            print(f"{ct}({len(coeff)-1}) -> ", end='')
    
def remez_recursion(routes: list, data: remezData, p_num: int, e_num: int, max_n: int, eb: EB, approx_mode: str, printmode="debug") -> list:
    n = max_n-1 if approx_mode == "odd" else max_n
    intervals = data.interval_log[-1]

    # 현재 경로가 완료한 경로보다 성능이 낮을경우 반환
    if len(data.pre_log) >= 2:
        pre_intervals = data.interval_log[-2]
        len_int = sum(i[1]-i[0] for i in intervals)
        len_int2 = sum(i[1]-i[0] for i in pre_intervals)
        if len_int2 < len_int:
            # print("interval length False")
            return routes
        
    if len(routes) > 0:
        if routes[0].compare(data, "onRoute") is False:
            # print("route complexity False")
            return routes
        
    if intervals[0][1] > intervals[1][0]:
            return routes
        
    def evalF(x):
        return 1 if x >= 0.5 else 0
    
    if approx_mode == "cl":
        coeff = [0, 0, 3, -2]
        max_err, next_intervals = calculate_next_remez(coeff, evalF, eb, intervals, 2)
        try:
            routes, _ = end_route(max_err, data, coeff, next_intervals, p_num, e_num, max_n, eb, routes, printmode)
        except Exception as e:
            print()
        
    else:
        while n > 0:
            coeff, max_err, next_intervals = remez_algorithm(n, intervals, evalF, approx_mode, eb, 2, "normal")
            try:
                routes, flag = end_route(max_err, data, coeff, next_intervals, p_num, e_num, max_n, eb, routes, printmode)
            except Exception as e:
                print()
            if flag is False:
                # break
                pass
            n -= 2
    return routes

def end_route(max_err: float, data: remezData, coeff: list, next_intervals: list, p_num: int, e_num: int, max_n: int, eb: EB, routes: list, print_mode: str):
    flag = False
    pre = -int(log2(max_err))
    if pre < 0 or pre < data.pre_log[-1]:
        return routes, flag

    # 새로운 remezData 객체 생성 후 정보 갱신
    newData = copy.deepcopy(data)
    newData.update(coeff, pre, next_intervals)

    if pre >= e_num: # 성공
        if len(routes) == 0:
            # print("Init!")
            routes = [newData]
            flag = True
            # write_params(p_num, newData)
        else:
            comp_result = routes[0].compare(newData)
            if comp_result is True:
                # print("Change!")
                routes = [newData]
                # newData.print_params()
                # write_params(p_num, newData)
                flag = True
            elif comp_result == 99:
                # print("Add!")
                routes.append(newData)
                # newData.print_params()
                flag = True
            # else:
            #     print("NO!")

        
        debug_print(f"n={len(coeff)-1}, pre {pre}, step {len(data.coeff_log) + 1}, intervals: ", print_mode, "")
        if print_mode == "debug":
            print_intervals(next_intervals, 10)
            
        return routes, flag
    else:
        routes = remez_recursion(routes, newData, p_num, e_num, max_n, eb, "even", print_mode)
        routes = remez_recursion(routes, newData, p_num, e_num, max_n, eb, "odd", print_mode)
        routes = remez_recursion(routes, newData, p_num, e_num, max_n, eb, "cl", print_mode)
    
    return routes, flag
    

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
            
        # 재귀 - 홀수/짝수/Cleanse 3가지 모두 진행
        routes = remez_recursion(routes, data, p_num, e_num, max_n, eb, "even", printmode)
        routes = remez_recursion(routes, data, p_num, e_num, max_n, eb, "odd", printmode)
        routes = remez_recursion(routes, data, p_num, e_num, max_n, eb, "cl", printmode)
    
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
            for coeff in routes[0].coeff_log:
                f.write(coeff2txt(coeff))
                f.write('\n')

    else: # 실패
        print("Approximation Failed.")
        # routes[0].print_params()
        
        
# def _mode_worker(args):
#     n, mode, intervals, p_num, e_num, max_n, eb, printmode = args

#     # evalF는 워커 내부에서 정의 (피클 이슈 회피를 위해 워커 외부 참조를 최소화)
#     def evalF(x):
#         return 1 if x < 0.5 else 0

#     # 초기 근사 (원 코드와 동일하게 "even"으로 초기화)
#     coeff, max_err, next_intervals = remez_algorithm(n, intervals, evalF, "even", eb, 2, "normal")
#     if max_err == 99:
#         # 실패 시 빈 리스트 반환
#         return []

#     pre = -int(log2(max_err))

#     # remezData 생성 및 갱신
#     data = remezData(intervals)
#     data.update(coeff, pre, next_intervals)

#     # 해당 (n, mode)에 대한 경로 리스트 반환
#     # routes 인자는 누적용이므로 개별 작업에선 빈 리스트로 시작
#     return remez_recursion([], data, p_num, e_num, max_n, eb, mode, printmode)


# def multi_remez_MT(p_num: int, e_num: int, max_n: int, eb, printmode: str):
#     p = pow(2, p_num)
#     e = eb.Bc / eb.scale
#     intervals = [[0, e], [(1 - e), (p - 1 + e)]]

#     modes = ["even", "odd", "cl"]
#     ns = list(range(2, max_n + 1, 2))
#     tasks = [(n, mode, intervals, p_num, e_num, max_n, eb, printmode) for n in ns for mode in modes]

#     # 워커(프로세스) 수를 32로 고정
#     with Pool(processes=32) as pool:
#         results = pool.map(_mode_worker, tasks, chunksize=1)

#     routes = []
#     for r in results:
#         routes.extend(r)

#     if len(routes) > 0 and routes[0].pre_log[-1] >= e_num:
#         print(f"Total {len(routes)} routes.")
#     compare_final(p_num, e_num, routes, eb)

# def compare_final(p_num, e_num, routes, eb):
#     # 최종 비교
#     bob = [routes[0]]

#     for i in range(1, len(routes)):
#         route = routes[i]
#         res = bob[0].compare(route, "end")
#         if res is True:
#             bob = [route]
#         elif res == 99:
#             bob.append(route)

#     bob[0].print_params()
    
#     # 텍스트파일 저장
#     with open(f"doc/coeff_{p_num}_{e_num}.txt", "w", encoding="utf-8") as f:
#         for coeff in bob[0].coeff_log:
#             f.write(coeff2txt(coeff))
#             f.write('\n')
    
#     # print(f"(p, e): {p_num}, {e_num}")
#     # # for b in bob:
#     #     # b.print_params()
#     # bob[-1].print_params()
    
    
# def write_params(p_num, newData, doc_dir="doc"):
#     """
#     기존 doc/ 위치에 coeff_{p_num}_{depth}.txt 파일들을 검사한 뒤,
#     - 현재 depth보다 작은 파일이 하나라도 있으면 아무 작업도 하지 않음.
#     - 현재 depth보다 큰 파일이 있으면 모두 삭제.
#     """
#     depth = newData.total_CalData.depth
#     file_pattern = re.compile(rf"^coeff_{p_num}_(\d+)\.txt$")

#     existing_files = [f for f in os.listdir(doc_dir) if file_pattern.match(f)]

#     # 기존 depth들 추출
#     existing_depths = []
#     for f in existing_files:
#         m = file_pattern.match(f)
#         if m:
#             existing_depths.append(int(m.group(1)))

#     # 더 작은 depth 파일이 있는지 확인
#     smaller_exists = any(d <= depth for d in existing_depths)
#     if smaller_exists:
#         return

#     # 현재 depth보다 큰 파일 삭제
#     for d in existing_depths:
#         if d > depth:
#             target = os.path.join(doc_dir, f"coeff_{p_num}_{d}.txt")
#             try:
#                 os.remove(target)
#                 # print(f"삭제됨: {target}")
#             except FileNotFoundError:
#                 pass

#     # 새 파일 생성
#     filepath = os.path.join(doc_dir, f"coeff_{p_num}_{depth}.txt")
#     with open(filepath, "w", encoding="utf-8") as f:
#         for coeff in newData.coeff_log:
#             f.write(coeff2txt(coeff))
#             f.write('\n')
            
# def _merge_routes(base_routes, new_routes):
#     """
#     base_routes와 new_routes를 remezData.compare("end") 규칙으로 병합
#     """
#     if not new_routes:
#         return base_routes
#     if not base_routes:
#         return new_routes

#     merged = base_routes[:]  # shallow copy
#     for route in new_routes:
#         res = merged[0].compare(route, "end")
#         if res is True:
#             merged = [route]
#         elif res == 99:
#             merged.append(route)
#         # res가 False면 무시
#     return merged

# # --- end_route 내부 병렬화 (재귀 분기 even/odd/cl 병렬 실행) ---
# def end_route_MT(max_err: float, data: remezData, coeff: list, next_intervals: list,
#               p_num: int, e_num: int, max_n: int, eb: EB, routes: list, print_mode: str):
#     flag = False
#     pre = -int(log2(max_err))
#     if pre < 0 or pre < data.pre_log[-1]:
#         return routes, flag

#     # 새로운 remezData 객체 생성 후 정보 갱신
#     newData = copy.deepcopy(data)
#     newData.update(coeff, pre, next_intervals)

#     if pre >= e_num: # 성공
#         if len(routes) == 0:
#             routes = [newData]
#             flag = True
#             write_params(p_num, newData)
#         else:
#             comp_result = routes[0].compare(newData)
#             if comp_result is True:
#                 routes = [newData]
#                 write_params(p_num, newData)
#                 flag = True
#             elif comp_result == 99:
#                 routes.append(newData)
#                 flag = True

#         debug_print(f"n={len(coeff)-1}, pre {pre}, step {len(data.coeff_log) + 1}, intervals: ", print_mode, "")
#         if print_mode == "debug":
#             print_intervals(next_intervals, 10)
#         return routes, flag

#     else:
#         # ---- 여기부터 내부 멀티스레딩: 세 가지 분기를 병렬로 탐색 ----
#         modes = ("even", "odd", "cl")
#         # 각 분기에서 공용 routes를 동시에 수정하지 않도록, 빈 routes로 시작해서 결과만 병합
#         partial_results = []

#         # ThreadPoolExecutor를 사용 (프로세스 내에서 안전, GIL-해제되는 연산/넘파이 비중 시 유효)
#         # 만약 remez_algorithm이 순수 파이썬 CPU 바운드라면, ThreadPool 대신 외부 Pool에서 확장 권장
#         with ThreadPoolExecutor(max_workers=len(modes)) as ex:
#             futures = [
#                 ex.submit(remez_recursion, [], newData, p_num, e_num, max_n, eb, mode, print_mode)
#                 for mode in modes
#             ]
#             for fut in as_completed(futures):
#                 try:
#                     partial = fut.result()
#                 except Exception:
#                     partial = []
#                 partial_results.append(partial)

#         # 병합
#         for pr in partial_results:
#             routes = _merge_routes(routes, pr)

#         return routes, flag
