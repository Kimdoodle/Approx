import numpy as np
from math import log2
from .print import debug_print
from scipy.optimize import brentq
from cal_depth import check_coeff_type

''' 기타 필요한 함수들 '''
# 다항함수 평가
def evalP(coeff, x):
    return sum(coeff[i] * x**i for i in range(len(coeff)))

def error_abs(coeff, x, y):
    return abs(evalP(coeff, x) - y)

def error_func(coeff, x, y):
    return evalP(coeff, x) - y

# 구간 [a, b]에서 chebyshev node를 사용하여 m개의 점 샘플링
def sample_points(a, b, m) -> np.ndarray:
    if m <= 0:
        print("Sample 0 points!!!")
        return np.ndarray([], dtype=float)
    k = np.arange(m)
    x = np.cos(np.pi * (2 * k + 1) / (2 * m))
    return 0.5 * (b + a) + 0.5 * (b - a) * x

def sample_points_multi(m, intervals, print_mode) -> np.ndarray:
    # 구간의 총 길이를 계산
    lengths = np.array([end - start for start, end in intervals])
    total_length = lengths.sum()
    n = len(intervals)
    m_alloc = np.ones(n, dtype=int)
    left = m - n

    if total_length == 0:
        # 총 길이가 0이면 m을 균등하게 분배
        base = m // n
        remainder = m % n
        m_alloc = np.full(n, base, dtype=int)
        if remainder > 0:
            m_alloc[:remainder] += 1
    else:
        # 각 구간의 길이 비율로 m을 분배, 남은 점은 긴 길이부터 추가로 할당
        m_float = left * lengths / total_length
        additional = np.floor(m_float).astype(int)
        m_alloc += additional

        remainder = m - m_alloc.sum()
        if remainder > 0:
            idx_desc = np.argsort(-lengths)
            for idx in idx_desc[:remainder]:
                m_alloc[idx] += 1

    # 구간 별 샘플링
    chunks = []
    for (s, e), num in zip(intervals, m_alloc):
        debug_print(f"sample {num} points in interval[{s}, {e}]", print_mode)
        if num > 0:
            chunks.append(sample_points(s, e, num))

    x = np.concatenate(chunks)
    x.sort()
    return x

def generate_points(start: np.float64, end: np.float64, step=np.float64(1e-10)):
    if start == end:
        return np.array([np.float64(start)], dtype=np.float64)

    max_points = 10000  # 최대 허용 포인트 수
    step = np.float64(step)

    while True:
        # 예상되는 포인트 개수 (양 끝 포함)
        size_est = int(np.floor((end - start) / step)) + 1

        # 포인트 개수가 허용 범위 내면 실제 생성
        if size_est <= max_points:
            # 직접 계산 (누적오차 방지)
            points = start + np.arange(size_est, dtype=np.float64) * step

            # 마지막 점이 end보다 작으면 end를 추가하여 양 끝 포함 보장
            if points[-1] < end:
                points = np.append(points, end)
            # 오차로 end를 살짝 초과한 경우 end로 클램프
            elif points[-1] > end:
                points[-1] = end

            return points

        # 초과 시 step 증가 (샘플링 간격 10배 확대)
        step *= np.float64(10.0)

# 구간 [a, b]에서 f(x)-p(x)의 근을 계산
def find_intersection(coeff, evalF, a, b, tol=1e-9) -> list:
    x_points = generate_points(np.float64(a), np.float64(b))
    cand = []

    def f(x):
        return evalP(coeff, x) - evalF(x)

    for i in range(len(x_points) - 1):
        x1, x2 = x_points[i], x_points[i+1]
        y1, y2 = f(x1), f(x2)

        # 거의 0인 경우(경계 중복 방지: 왼쪽 끝만 채택)
        if np.isfinite(y1) and np.isclose(y1, 0.0, atol=tol):
            cand.append(float(x1))
            continue
        if np.isfinite(y2) and np.isclose(y2, 0.0, atol=tol):
            # x2는 다음 구간의 x1이므로 여기서는 추가하지 않음
            pass
        elif np.isfinite(y1) and np.isfinite(y2) and y1 * y2 < 0:
            try:
                root = brentq(f, x1, x2)
                cand.append(float(root))
            except ValueError:
                # 브래킷 문제나 불연속 등은 스킵
                pass

    # 허용오차로 중복 제거
    cand.sort()
    roots = []
    for x in cand:
        if not roots or not np.isclose(x, roots[-1], atol=tol):
            roots.append(x)
    return roots

# 구간 변경
def slice_interval(approx_mode: str, intervals: list) -> list:
    # if approx_mode in ["all", "etc"]:
    #     return intervals
    
    # 기존 구간 intervals 정렬
    new_intervals = []
    for start, end in intervals:
        if start < 0 and end < 0:
            continue
        elif start < 0 and end > 0:
            if abs(start) > abs(end):
                new_intervals.append([0.0, abs(start)])
            else:
                new_intervals.append([0.0, end])
        else:
            new_intervals.append([start, end])

    return new_intervals
            
# 행렬식 구성
def create_matrix(powers: list, x: np.ndarray, evalF) -> tuple[np.ndarray, np.ndarray]:
    A_matrix = []
    y_matrix = []
    for xindex in range(len(x)):
        x_sample = x[xindex]
        rowA = [x_sample ** p for p in powers] + [(-1) ** xindex]
        y = evalF(x_sample)
        A_matrix.append(rowA)
        y_matrix.append(y)

    A_matrix = np.array(A_matrix, dtype=np.float64)
    y_matrix = np.array(y_matrix, dtype=np.float64)
    return A_matrix, y_matrix


# 행렬식 계산 및 계수벡터 연산
def solve_matrix(A: np.ndarray, y: np.ndarray, n: int, powers: list) -> tuple[list, float]:
    try:
        B = np.linalg.solve(A, y)
    except np.linalg.LinAlgError:
        # print("singular matrix error!")
        # print(f"A:\n{A}")
        return [-1], -1
    E = B[-1]
    coeff = []
    for k in range(n + 1):
        if k in powers:
            coeff.append(B[0])
            B = B[1:]
        else:
            coeff.append(0.0)
    return coeff, E

# 지역 극점의 x좌표 계산
def calculate_local_max(coeff, evalF, intervals: list) -> tuple:
    cross_interval = [] # (start, end) tuple로 구성된 list
    max_point_x = []
    max_point_y = []
    
    # 근을 기준으로 새로운 구간을 정의
    for start, end in intervals:
        roots = find_intersection(coeff, evalF, start, end)
        if len(roots) == 0:
            cross_interval.append((np.float64(start), np.float64(end)))
        else:
            if start not in roots:
                roots.insert(0, start)
            if end not in roots:
                roots.append(end)
            for i in range(len(roots)-1):
                cross_interval.append((roots[i], roots[i+1]))

    for start, end in cross_interval:
        try:
            # local_x = np.linspace(start, end, 1000)
            local_x = generate_points(start, end)
            local_y = [error_abs(coeff, x, evalF(x)) for x in local_x]
            max_index = np.argmax(local_y)
            max_x = local_x[max_index]
            max_y = local_y[max_index]
            max_point_x.append(max_x)
            max_point_y.append(max_y)
        except Exception as e:
            pass
    
    return max_point_x, max_point_y

# 종료조건 판단
def decide_exit(errors: list, threshold: float, print_mode: str) -> bool:
    max_err, min_err = max(errors), min(errors)
    if abs(max_err - min_err) <= 1e-9:
        return True
    debug_print(f"Min/Max error\t\t: {min_err:.5f}, {max_err:.5f}", print_mode)
    if min_err != 0:
        exit_err = (max_err - min_err) / min_err
    else:
        exit_err = 0.0
    debug_print(f"Exit error\t\t: {exit_err}", print_mode)

    if exit_err <= threshold:
        return True
    else:
        return False

# remez 다음을 위한 interval, max_err 계산
# mode1 : errbound계산에 B_clean 고려
# mode2 : errbound계산에 B_scale만 고려
def calculate_next_remez(coeff, evalF, eb, intervals):
    # if (
    #     len(coeff) >= 3
    #     and coeff[0] != 0
    #     and coeff[-1] != 0
    #     and coeff[-2] != 0
    #     and all(x == 0 for i, x in enumerate(coeff) if i not in (0, len(coeff)-1, len(coeff)-2))
    # ):
    #     pass
    ni = []
    max_err0 = 0
    max_err1 = 0
    coeff_type = check_coeff_type(coeff)
    
    for i, interval in enumerate(intervals):
        start = interval[0]
        end = interval[1]
        if start == end:
            val = evalF(start) - evalP(coeff, start)
            err2 = eb.cal_bound(start, coeff, coeff_type)/eb.scale
            val1 = val + err2
            val2 = val - err2
            ni.append([val2, val1])
            if i == 0:
                max_err0 = max(val1, max_err0)
            else:
                max_err1 = max(val1, max_err1)
        else:
            points = generate_points(start, end)
            try:
                p_vals1      = np.fromiter((evalP(coeff, p) + eb.cal_bound(p, coeff, coeff_type)/eb.scale for p in points), dtype=float, count=points.size)
                p_vals2      = np.fromiter((evalP(coeff, p) - eb.cal_bound(p, coeff, coeff_type)/eb.scale for p in points), dtype=float, count=points.size)
                f_vals      = np.fromiter((evalF(p) for p in points), dtype=float, count=points.size)

                err_vals = np.concatenate((np.abs(f_vals - p_vals1), np.abs(f_vals - p_vals2)))
                vals = np.concatenate((p_vals1, p_vals2))
                ni.append([np.min(vals), np.max(vals)])
                
                # 최대오차
                if i == 0:
                    max_err0 = max(np.max(err_vals), max_err0)
                else:
                    max_err1 = max(np.max(err_vals), max_err1)
            except Exception as e:
                print("CNR err")
                print(coeff)
                print(e)
                return 99, 99, [[0,0], [1,1]]
    
    # 정렬 - 구간이 긴 쪽이 다음 근사에서 0으로 수렴해야함
    # if (ni[0][1]-ni[0][0]) < (ni[1][1]-ni[1][0]):
    #     ni.reverse()

    if ni[0][0] > ni[1][0]:
        ni.reverse()

    return max_err0, max_err1, ni




