from src_remez.algorithm import remez_algorithm
from src_remez.print_plot import *
from src_remez.print import print_intervals, debug_print, coeff2txt
from src_approx.approx_remez import remezData, multi_remez
from src_remez.basic import calculate_next_remez
from src_errbound.error_bound import EB
# from parse import parse_into_remezData, get_file
from src_approx.util_MT import compare_depth_route, _STORE, _STORE_LOCK
from math import log2
import os, copy

from src_approx.util_MT import (
    _ensure_exec, _inc, _dec, wait_for_all, shutdown_exec,
    shared_init, shared_get, shared_update_if_better
)
from src_approx.file_MT import init_lockfile, read_params_threadsafe, write_params_threadsafe, write_params

# ====== CPU affinity(0~31) ======
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
try:
    os.sched_setaffinity(0, set(range(32)))
except Exception:
    pass

# Remez - Step1
def multi_remez(p_num: int, e_num: int, max_n: int, eb: EB, printmode: str, multiThread: bool):

    p = pow(2, p_num)
    e = eb.Bc / eb.scale
    intervals = [[0, e], [(1 - e), (p - 1 + e)]]

    def evalF(x):
        return 1 if x < 0.5 else 0

    for n in range(2, max_n + 1, 2):
        try:
            coeff, max_err0, max_err1, next_intervals = remez_algorithm(n, intervals, evalF, "even", eb, 2, "normal")
        except Exception:
            continue
        if max_err0 == 99:
            debug_print(f"Init approximation n={n} failed.", printmode)
            continue

        # pre0, pre1 = -log2(max_err0), -log2(max_err1)
        debug_print(f"Init approximation\nn={n}, step 0, intervals: ", printmode, "")
        if printmode == "debug":
            print_intervals(next_intervals, 10)

        data = remezData(intervals)
        # data.update(coeff, (pre0, pre1), next_intervals)

        end_route(data, coeff, max_err0, max_err1, next_intervals, p_num, e_num, max_n, eb, printmode, multiThread)
        # remez_recursion(data, p_num, e_num, max_n, eb, "even", printmode, multiThread)
        # remez_recursion(data, p_num, e_num, max_n, eb, "odd",  printmode, multiThread)
        # remez_recursion(data, p_num, e_num, max_n, eb, "cl", printmode, multiThread)
        # remez_recursion(data, p_num, e_num, max_n, eb, "all", printmode, multiThread)
        

        # 모든 프로세스 마무리
        wait_for_all()
        
        # 데이터 기록 
        with _STORE_LOCK:
            for depth in sorted(_STORE.keys()):
                data = _STORE[depth]
                pre = data.pre_log[-1]
                if pre[0] >= e_num and pre[1] >= e_num: # 성공
                    doc_dir = "doc/"
                    filepath = os.path.join(doc_dir, f"coeff_{p_num}_{e_num}.txt")
                    with open(filepath, "w", encoding="utf-8") as f:
                        for coeff in data.coeff_log:
                            f.write(coeff2txt(coeff))
                            f.write('\n')
                    break
        # best_data = shared_get(p_num)
        # if best_data is not None:
        #     write_params(p_num, e_num, best_data)
        #     debug_print(f"[write_params] Final best remezData saved for p_num={p_num}", printmode)
        # else:
        #     debug_print(f"[write_params] No best data found for p_num={p_num}", printmode)

        shutdown_exec()

            

# Step 2~
def remez_recursion(data: remezData, p_num: int, e_num: int, max_n: int, eb: EB, approx_mode: str, printmode: str, multiThread: bool):
    try:
        n = max_n - 1 if approx_mode == "odd" else max_n
        intervals = data.interval_log[-1]

        if approx_mode == "cl":  
            try:
                step2_CL(intervals, eb, data, p_num, e_num, max_n, printmode, multiThread)
            except Exception as e:
                print("Step2_CL error")
                print(e)
                return
        else:
            try:
                step2_remez(intervals, approx_mode, eb, n, data, p_num, e_num, max_n, printmode, multiThread)
            except Exception as e:
                print(f"Step2_Remez {approx_mode} error")
                print(e)
                return
            
    except Exception as e:
        print("Recursion error")
        print(e)
        return

def step2_CL(intervals, eb: EB, newData: remezData, p_num: int, e_num: int, max_n: int, print_mode: str, multiThread: bool):
    if multiThread:
        for cur_n in range(7, 2, -2):
            _inc()
            fut = _EXEC_submit(cal_CL, intervals, eb, cur_n)
            fut.add_done_callback(lambda f, d=newData: (_dec(), end_route(d, *f.result(), p_num, e_num, max_n, eb, print_mode, multiThread)))
    else:
        for cur_n in range(7, 2, -2):
            coeff, max_err0, max_err1, next_interval = cal_CL(intervals, eb, cur_n)
            if coeff[0] == 1.0 and coeff[-1] == 0.0:
                pass
            end_route(newData, coeff, max_err0, max_err1, next_interval, p_num, e_num, max_n, eb, print_mode, multiThread)


def step2_remez(intervals, approx_mode: str, eb: EB, n_start: int, newData: remezData, p_num: int, e_num: int, max_n: int, print_mode: str, multiThread: bool):
    minus = -2 if approx_mode in ["odd", "even"] else -1
    if multiThread:
        _ensure_exec()
        for cur_n in range(n_start, 1, minus):
            if cur_n == 1:
                print("!!!!!!!!!!!!!!!!!!!!!!!!")
            _inc()
            fut = _EXEC_submit(cal_remez, cur_n, intervals, approx_mode, eb)
            fut.add_done_callback(lambda f, d=newData: (_dec(), end_route(d, *f.result(), p_num, e_num, max_n, eb, print_mode, multiThread)))
    else:
        for cur_n in range(n_start, 1, minus):
            coeff, max_err0, max_err1, next_interval = cal_remez(cur_n, intervals, approx_mode, eb)
            if coeff[0] == 1.0 and coeff[-1] == 0.0:
                pass
            end_route(newData, coeff, max_err0, max_err1, next_interval, p_num, e_num, max_n, eb, print_mode, multiThread)


# 내부 제출 헬퍼 (util_MT의 실행자에 의존)
def _EXEC_submit(fn, *args, **kwargs):
    # util_MT에서 관리되는 _EXEC에 접근하기 위해 import 대신 지연 참조
    # (ProcessPoolExecutor 인스턴스는 util_MT._ensure_exec로 생성됨)
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
            case 7:
                coeff = [0, 0, 0, 0, 0, 0, 7, -6]
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

# 수렴결과 확인
def end_route(data: remezData, coeff, max_err0, max_err1, next_intervals, p_num: int, e_num: int, max_n: int, eb: EB, print_mode: str, multiThread: bool):
    try:
        pre0, pre1 = -log2(max_err0), -log2(max_err1)
        
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
            4. (test) 현재 사용한 depth의 최적 결과보다 두 정밀도가 낮은 경우
        '''
        if compare_depth_route(data) is False:
            return

        # if pre0 >= e_num and pre1 >= e_num:
        #     # write_params_threadsafe(p_num, newData)
        #     shared_update_if_better(p_num, newData)
        #     debug_print(f"n={len(coeff)-1}, step {len(data.coeff_log)+1}, depth: {newData.total_CalData.depth}", print_mode, '\n')
        #     return

        # for m in ["odd", "even", "cl"]:
        m = "all"
        remez_recursion(newData, p_num, e_num, max_n, eb, m, print_mode, multiThread)

    except Exception as e:
        print("endroute err")
        print(e)
        return
 