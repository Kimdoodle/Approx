import numpy as np
from .basic import *
from .basic_multi import find_best_combination, find_best_combination_v2, find_best_combination_v3
from .print import *
from .print_plot import *
from typing import Literal
import sympy as sp
from numpy.polynomial import polynomial as P
from src_errbound.error_bound import EB

# remez1 - errbound 계산에 초기 B_clean을 고려함
# remez2 - errbound 계산에 B_scale만을 고려함
def remez_algorithm(n: int, intervals: list, evalF, approx_mode: str, eb: EB, err_mode: int, print_mode: str) -> tuple:  
    intervals = slice_interval(approx_mode, intervals)

    # 1: x 샘플링
    if approx_mode == "odd":
        sample_number = int((n + 1) / 2) + 1
    elif approx_mode == "even":
        sample_number = int(n / 2 + 2)
    elif approx_mode == "all":
        sample_number = n + 2
    x_samples = sample_points_multi(sample_number, intervals, print_mode)

    for _ in range(50):
        # 2: 행렬식 구성
        if approx_mode == "even":
            powers = [i for i in range(n % 2, n + 1, 2)]
        elif approx_mode == "all":
            powers = [i for i in range(n+1)]
        A_matrix, y_matrix = create_matrix(powers, x_samples, evalF)

        # 3: 행렬식 연산
        try:
            assert A_matrix.shape[0] == A_matrix.shape[1], \
                    f"A가 직사각 행렬입니다. shape={A_matrix.shape[0]}, {A_matrix.shape[1]}"
        except Exception as e:
            # print("행렬식 연산 오류")
            # print(e)
            return [-1], 99, []
        try:
            coeff, E = solve_matrix(A_matrix, y_matrix, n, powers)
        except Exception as e:
            # print("??")
            # print(e)
            return [-1], 99, []

        # 4: 지역 극값 연산
        max_point_x, max_point_y = calculate_local_max(coeff, evalF, intervals)

        if len(max_point_x) == sample_number:
            best_points_x = max_point_x
            best_points_y = max_point_y

        # 극점 개수가 부족한 경우
        elif len(max_point_x) < sample_number:
            req_points = sample_number - len(max_point_x)
            # print(f"Error\t\t\t: {req_points} more local max points required.")
            debug_print(f"Error\t\t\t: {req_points} more local max points required.", print_mode)
            # 각 구간의 양 끝점 중 오차가 큰 순서대로 삽입.
            boundary_data = []
            for start, end in intervals:
                boundary_data.append((start, error_abs(coeff, start, evalF(start))))
                boundary_data.append((end, error_abs(coeff, end, evalF(end))))
            boundary_data.sort(key= lambda x: x[1])

            added = 0
            for x_val, err in reversed(boundary_data):
                if added == req_points:
                    break
                x_np = np.float64(x_val)
                if x_np not in max_point_x:
                    max_point_x.append(x_np)
                    max_point_y.append(err)
                    added += 1

            # 2개의 점으로도 부족할 경우 오류를 출력하고 더미데이터 반환.
            if added < req_points:
                debug_print(f"Fatal Error\t\t: {req_points-added} more points need to be added.", print_mode)
                best_points_x = max_point_x
                best_points_y = max_point_y
                return [-1], 99, []
            else:
                max_point_x.sort()
                max_point_y.sort()
                best_points_x = max_point_x
                best_points_y = max_point_y

        # 극점 개수가 너무 많은 경우
        elif len(max_point_x) > sample_number:
            try:
                extra_points = len(max_point_x) - sample_number
                debug_print(f"{extra_points} more local max points exists.", print_mode)
                # best_points_x, best_points_y  = find_best_combination_v3(max_point_x, max_point_y, extra_points, coeff, evalF)
                best_points_x, best_points_y = find_best_combination(max_point_x, max_point_y, sample_number, coeff, evalF, print_mode)
            except TypeError as e:
                # print("극점 개수 오류")
                # print(e)
                return [-1], 99, []

        best_error_abs = [abs(x) for x in best_points_y]

        # 5: 종료조건 판단
        if decide_exit(best_error_abs, 1e-2, print_mode):
            try:
                max_err, next_intervals = calculate_next_remez(coeff, evalF, eb, intervals, err_mode)
                if next_intervals[0][1] < next_intervals[1][0]:
                    return coeff, max_err, next_intervals
                else:
                    debug_print("구간 대소관계 오류", print_mode)
                    return [-1], 99, []
            
            except Exception as e:
                print("구간 계산 오류")
                print(intervals)
                print(e)
                return [-1], 99, []

        else:
            x_samples = best_points_x

        debug_print("", print_mode)
    
    debug_print("Approx failed.", print_mode)
    return [-1], 99, []


def cleanse(eb: EB, intervals: list[list[float]]):
    def evalF(x):
        return 1 if x >= 0.5 else 0
    
    def evalP(x):
        return -2 * pow(x, 3) + 3 * pow(x, 2)
    
    ni = []
    max_err = 0.0
    for start, end in intervals:
        points = generate_points(np.float64(start), np.float64(end), np.float64(1e-8))
        p_vals1 = np.fromiter((evalP(p) + eb.cal_bound_cleanse(p)/eb.scale for p in points), dtype=float, count=len(points))
        p_vals2 = np.fromiter((evalP(p) - eb.cal_bound_cleanse(p)/eb.scale for p in points), dtype=float, count=len(points))
        f_vals = np.fromiter((evalF(p) for p in points), dtype=float, count=len(points))
        
        err_vals = np.concatenate((np.abs(f_vals - p_vals1), np.abs(f_vals - p_vals2)))
        vals = np.concatenate((p_vals1, p_vals2))
        ni.append([np.min(vals), np.max(vals)])
    
    # 최대오차
    max_err = max(np.max(err_vals), max_err)
    
    if (ni[0][1]-ni[0][0]) < (ni[1][1]-ni[1][0]):
        ni.reverse()

    return max_err, ni
    
            