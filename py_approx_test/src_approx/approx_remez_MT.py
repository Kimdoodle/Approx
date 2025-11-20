from src_remez.algorithm import remez_algorithm
from src_remez.print_plot import *
from src_remez.print import print_intervals, debug_print, write_params
from src_approx.approx_remez import remezData
from src_remez.basic import calculate_next_remez
from cal_depth import check_coeff_type
from src_errbound.error_bound import EB
from math import log2
import os, copy, time

from src_approx.util_MT import (
    _EXEC, shared_clear, _ensure_manager, _get_pre_lock,
    wait_for_all, shutdown_exec, _inc, _dec, _ensure_exec,
)
import src_approx.util_MT as util_MT

# ====== CPU affinity(0~31) ======
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
try:
    os.sched_setaffinity(0, set(range(32)))
except Exception:
    pass

# Remez - Step1
def multi_remez(p_num: int, e_num: int, max_n: int, eb: EB, tail: str, printmode: str, multiThread: bool, useDP: bool):
    shared_clear()
    p = pow(2, p_num)
    e = eb.Bc / eb.scale
    intervals = [[0, e], [(1 - e), (p - 1 + e)]]

    def evalF(x):
        return 1 if x < 0.5 else 0

    for n in range(max_n, 1, -2):
        try:
            coeff, max_err0, max_err1, next_intervals = remez_algorithm(n, intervals, evalF, "even", eb, 2, "normal")
        except Exception:
            continue
        if max_err0 == 99 and max_err1 == 99:
            debug_print(f"Init approximation n={n} failed.", printmode)
            continue

        debug_print(f"Init approximation\nn={n}, step 0, intervals: ", printmode, "")
        if printmode == "debug":
            print_intervals(next_intervals, 10)

        data = remezData(intervals)
        end_route(data, coeff, max_err0, max_err1, next_intervals, p_num, e_num, max_n, eb, tail, printmode, multiThread, useDP)

        # 모든 프로세스 마무리
        wait_for_all()
        shutdown_exec()

# 수렴결과
def end_route(data: remezData, coeff, max_err0, max_err1, next_intervals, p_num: int, e_num: int, max_n: int, eb: EB, tail: str, print_mode: str, multiThread: bool, useDP: bool):
    try:
        pre0 = safe_neg_log2(max_err0)
        pre1 = safe_neg_log2(max_err1)

        # 경로 지속 유효성 확인
        '''
            1. 두 정밀도가 모두 음수일 경우
            2. 두 정밀도가 모두 직전 정밀도보다 낮을 경우
            3. 두 구간이 겹치는 경우
        '''
        if pre0 < 0 or pre1 < 0:
            return
        if data.pre_log and pre0 < data.pre_log[-1][0] and pre1 < data.pre_log[-1][1]:
            return
        if next_intervals[0][1] > next_intervals[1][0]:
            return
        newData = copy.deepcopy(data)
        newData.update(coeff, (pre0, pre1), next_intervals)

        '''
            4. CL로 수렴하지 못한 경우
            5. 현재 사용한 depth의 최적 결과보다 두 정밀도가 낮은 경우
               혹은 정밀도 조건을 만족하는 경우(탐색 종료)
        '''        
        if tail == "oe" and len(newData.coeff_log) >= 1:
            t1 = check_coeff_type(newData.coeff_log[-1])
            if not(pre0 >= e_num and pre1 >= e_num) and t1 == "cl":
                return
        
        # 교체
        pre = int(min(pre0, pre1, 30))        
        util_MT._ensure_manager()
        lock = util_MT._get_pre_lock(pre)
        with lock:
            prev = util_MT._STORE.get(pre, None)
        
        if prev is None:
            change_route(p_num, e_num, pre, newData, tail)
        else:
            complexity_res = prev.compare(newData)
            pre_th = pre0 >= e_num and pre1 >= e_num
            if useDP: # DP
                pp0, pp1 = prev.pre_log[-1]
                np0, np1 = newData.pre_log[-1]
                pre_res = pp0 > np0 and pp1 > np1
                if complexity_res is False and (pre_th is True or pre_res is True):
                    if print_mode == "debug" and pre_th is True:
                        print("DP false!!")
                        print(f"{pp0}, {pp1}, {np0}, {np1}")
                        print(f"depth\t{prev.total_CalData.depth} < {newData.total_CalData.depth}")
                        print(f"cmult\t{prev.total_CalData.cmult} < {newData.total_CalData.cmult}")
                    return
                
            if pre_th and complexity_res is True:
                change_route(p_num, e_num, pre, newData, tail) 
                return
        
    except Exception as e:
        print("endroute err")
        print(e)
        return
    
    if tail == "oe":
        for m in ["odd", "even", "cl"]:
            remez_recursion(newData, p_num, e_num, max_n, eb, m, tail, print_mode, multiThread, useDP)
    elif tail == "all":
        remez_recursion(newData, p_num, e_num, max_n, eb, "all", tail, print_mode, multiThread, useDP)

# Step 2~
def remez_recursion(data: remezData, p_num: int, e_num: int, max_n: int, eb: EB, approx_mode: str, tail: str, printmode: str, multiThread: bool, useDP: bool):
    try:
        n = max_n - 1 if approx_mode == "odd" else max_n
        intervals = data.interval_log[-1]

        if approx_mode == "cl":  
            try:
                step2_CL(intervals, eb, data, p_num, e_num, max_n, tail, printmode, multiThread, useDP)
            except Exception as e:
                print("Step2_CL error")
                print(e)
                return
        else:
            try:
                step2_remez(intervals, approx_mode, eb, n, data, p_num, e_num, max_n, tail, printmode, multiThread, useDP)
            except Exception as e:
                print(f"Step2_Remez {approx_mode} error")
                print(e)
                return
            
    except Exception as e:
        print("Recursion error")
        print(e)
        return

def step2_CL(intervals, eb: EB, newData: remezData, p_num: int, e_num: int, max_n: int, tail: str, print_mode: str, multiThread: bool, useDP: bool):
    start = 3
    end = 3
    
    # if intervals[0][1] < (1/6)**(1/5):
    #     end = 8
    # el
    if intervals[0][1] < (1/4)**(1/3):
        end = 6
    elif intervals[0][1] < 0.5:
        end = 4
            
    if multiThread:
        _ensure_exec()

        def _done_callback(f, d=newData):
            try:
                coeff, max_err0, max_err1, next_interval = f.result()
                end_route(d, coeff, max_err0, max_err1, next_interval, p_num, e_num, max_n, eb, tail, print_mode, multiThread, useDP)
            finally:
                _dec()
            
        for cur_n in range(start, end, 2):
            _inc()
            fut = _EXEC_submit(cal_CL, intervals, eb, cur_n)
            fut.add_done_callback(_done_callback)
    else:
        for cur_n in range(start, end, 2):
            coeff, max_err0, max_err1, next_interval = cal_CL(intervals, eb, cur_n)
            if coeff[0] == 1.0 and coeff[-1] == 0.0:
                pass
            end_route(newData, coeff, max_err0, max_err1, next_interval, p_num, e_num, max_n, eb, tail, print_mode, multiThread, useDP)

def step2_remez(intervals, approx_mode: str, eb: EB, n_start: int, newData: remezData, p_num: int, e_num: int, max_n: int, tail: str, print_mode: str, multiThread: bool, useDP: bool):
    minus = -2 if approx_mode in ["odd", "even"] else -1
    if multiThread:
        _ensure_exec()

        def _done_callback(f, d=newData):
            try:
                coeff, max_err0, max_err1, next_interval = f.result()
                end_route(d, coeff, max_err0, max_err1, next_interval, p_num, e_num, max_n, eb, tail, print_mode, multiThread, useDP)
            finally:
                _dec()

        for cur_n in range(n_start, 1, minus):
            _inc()
            fut = _EXEC_submit(cal_remez, cur_n, intervals, approx_mode, eb)
            fut.add_done_callback(_done_callback)
    else:
        for cur_n in range(n_start, 1, minus):
            coeff, max_err0, max_err1, next_interval = cal_remez(cur_n, intervals, approx_mode, eb)
            if coeff[0] == 1.0 and coeff[-1] == 0.0:
                pass
            end_route(newData, coeff, max_err0, max_err1, next_interval, p_num, e_num, max_n, eb, tail, print_mode, multiThread, useDP)


# 내부 제출 헬퍼
def _EXEC_submit(fn, *args, **kwargs):
    from src_approx.util_MT import _EXEC  # type: ignore
    return _EXEC.submit(fn, *args, **kwargs)

# 계산
def cal_CL(intervals, eb: EB, n: int) -> tuple[list, float, float, list]:
    try:
        def evalF(x):
            return 1 if x >= 0.5 else 0

        match n:
            case 3:
                coeff = [0, 0, 3, -2]
            case 5:
                coeff = [0, 0, 0, 0, 5, -4]
            # case 7:
            #     coeff = [0, 0, 0, 0, 0, 0, 7, -6]
        if intervals == [[np.float64(-4.468064867817073e-22), np.float64(5.3653120512813356e-08)], [np.float64(0.9999999990211308), np.float64(0.9999999990211308)]]:
            print("")
        max_err0, max_err1, next_intervals = calculate_next_remez(coeff, evalF, eb, intervals)
        return coeff, float(max_err0), float(max_err1), next_intervals
    except Exception as e:
        print("cal_CL err")
        print(e)
        return [-1], 99, 99, []

def cal_remez(cur_n: int, intervals, approx_mode: str, eb: EB) -> tuple[list, float, float, list]:
    try:
        def evalF(x):
            return 1 if x >= 0.5 else 0
        coeff, max_err0, max_err1, next_intervals = remez_algorithm(cur_n, intervals, evalF, approx_mode, eb, 2, "normal")
        return coeff, float(max_err0), float(max_err1), next_intervals
    except Exception as e:
        print("cal_remez err")
        print(e)
        return [-1], 99, 99, []

def safe_neg_log2(x, eps=1e-12):
    x = max(x, eps)
    return -log2(x)

def compare_route(p_num, e_num, prev: remezData, newData: remezData, tail: str):
    # 동일한 수준의 데이터와 비교 시 계산복잡도, 정밀도가 모두 낮을 때 탐색을 종료(False반환)
    try:
        '''
            계산복잡도
            True: newData가 더 좋음
            False: 기존(prev)이 더 좋음
            99: 두 계산복잡도가 완전히 동일함
        '''
        complexity_res = prev.compare(newData)
        
        '''
            정밀도
            True: newData가 더 높음
            False: 기존이 더 높음
            99: 둘 중 하나씩 더 높음
        '''
        pp0, pp1 = prev.pre_log[-1]
        np0, np1 = newData.pre_log[-1]
        
        if pp0 < np0 and pp1 < np1:
            pre_res = True
        elif pp0 > np0 and pp1 > np1:
            pre_res = False
        else:
            pre_res = 99

        if complexity_res is False and pre_res is False:
            return False

        return True
    except Exception as e:
        print("compare route err")
        print(e)
        return
    
def change_route(p_num: int, e_num: int, pre: int, newData: remezData, tail: str):
    util_MT._ensure_manager()
    lock = util_MT._get_pre_lock(pre)
    
    try:
        with lock:
            util_MT._STORE[pre] = newData
            if pre == e_num:
                write_params(p_num, e_num, newData, tail)
                cd, cc, cp, ca = newData.total_CalData.depth, newData.total_CalData.cmult, newData.total_CalData.pmult, newData.total_CalData.cadd
                print(f"depth, cmult, pmult, add: {cd}, {cc}, {cp}, {ca}")
    
    except Exception as e:
        print("change route err")
        print(e)
        return