import numpy as np
from math import log2
from .print import debug_print

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

def generate_points(start: np.float64, end: np.float64, step=np.float64(1e-8)):
    max_points = 10000  # 최대 허용 포인트 수
    step = np.float64(step)

    while True:
        p = np.arange(start, end + step, step, dtype=np.float64)
        if p.size <= max_points:
            return p
        step *= np.float64(10.0)

# def sample_points_multi(m, intervals, print_mode) -> np.ndarray:
#     """
#     구간 비율에 따라 총 m개의 점을 여러 구간에 분배하여 샘플링합니다.
#     - 각 구간마다 최소 1개 이상 샘플링하도록 설계되어 있으므로 m >= len(intervals) 이어야 합니다.
#     - 모든 구간은 (start < end) 를 만족해야 합니다.
#     - intervals 는 비어있을 수 없습니다.
#     """
#     # === 입력/전제 조건 검증 ===
#     # m 타입/값 검증
#     if not isinstance(m, (int, np.integer)):
#         raise TypeError(f"'m' must be an integer, got {type(m).__name__}")
#     if m < 0:
#         raise ValueError(f"'m' must be >= 0, got {m}")

#     # intervals 검증 (비어있음, 구조/타입, start<end)
#     if not intervals:
#         raise ValueError("'intervals' must not be empty.")

#     try:
#         lengths = np.array([float(end) - float(start) for start, end in intervals], dtype=float)
#     except Exception as e:
#         raise TypeError("'intervals' must be an iterable of (start, end) numeric pairs.") from e

#     if lengths.ndim != 1 or lengths.size != len(intervals):
#         raise ValueError("'intervals' must be a flat list of (start, end) pairs.")

#     # 각 구간 길이 유효성: 반드시 양수
#     if np.any(~np.isfinite(lengths)):
#         raise ValueError("Interval bounds must be finite numbers.")
#     if np.any(lengths <= 0):
#         bad_idxs = np.where(lengths <= 0)[0].tolist()
#         raise ValueError(f"All intervals must satisfy end > start. Invalid indices: {bad_idxs},\n intervals: {intervals}")

#     # 총 길이 검증
#     total_length = float(lengths.sum())
#     if not np.isfinite(total_length) or total_length <= 0:
#         raise ValueError("Sum of interval lengths must be positive.")

#     n = len(intervals)

#     # 분배 설계상 각 구간에 최소 1개 배정 -> m은 구간 수 이상이어야 함
#     if m < n:
#         raise ValueError(f"'m' must be at least the number of intervals ({n}), got {m}")

#     # === 분배 계산 ===
#     # 각 구간에 최소 1개 할당 후 남은 개수를 길이 비율로 배분
#     m_alloc = np.ones(n, dtype=int)
#     left = m - n  # 남은 개수 (>= 0 보장)

#     # 비율 배분 (부동소수점 -> 내림)
#     m_float = left * (lengths / total_length)  # total_length > 0 보장
#     if not np.all(np.isfinite(m_float)):
#         raise RuntimeError("Internal error: non-finite allocation encountered.")

#     additional = np.floor(m_float).astype(int)
#     m_alloc += additional

#     # 내림으로 인해 남는 잔여 개수 재분배 (긴 구간부터)
#     remainder = int(m - int(m_alloc.sum()))
#     if remainder > 0:
#         idx_desc = np.argsort(-lengths)  # 길이 내림차순 인덱스
#         for idx in idx_desc[:remainder]:
#             m_alloc[idx] += 1

#     # === 구간별 샘플링 ===
#     chunks = []
#     for (s, e), num in zip(intervals, m_alloc):
#         debug_print(f"sample {num} points in interval[{s}, {e}]", print_mode)
#         if num > 0:
#             # sample_points(s, e, num) 은 외부에 정의되어 있다고 가정
#             chunks.append(sample_points(float(s), float(e), int(num)))

#     # 모든 배정이 0인 경우는 위의 검증/분배상 발생하지 않지만, 방어적 처리
#     if not chunks:
#         return np.array([], dtype=float)

#     x = np.concatenate(chunks)
#     x.sort()
#     return x

# 구간 [a, b]에서 f(x)-p(x)의 근을 계산
def find_intersection(coeff, evalF, a, b) -> list:
    # num_points = int((b - a) * 10000)
    # x_points = np.linspace(a, b, num_points)
    x_points = generate_points(a, b)
    p_vals = np.fromiter((evalP(coeff, x) for x in x_points), dtype=np.float64)
    f_vals = np.fromiter((evalF(x) for x in x_points), dtype=np.float64)
    # p_vals = np.array([evalP(coeff, x) for x in x_points])
    # f_vals = np.array([evalF(x) for x in x_points])
    errors = p_vals - f_vals
    signs = np.sign(errors)
    x_cross = []

    # 4) 정확히 0인 지점
    zero_idx = np.nonzero(errors == np.float64(0.0))[0]
    if zero_idx.size > 0:
        x_cross.extend(x_points[zero_idx].tolist())

    # 5) 인접한 포인트 간 부호 변화 구간
    change_idx = np.nonzero(signs[:-1] * signs[1:] < 0)[0]
    if change_idx.size > 0:
        x0 = x_points[change_idx]
        x1 = x_points[change_idx + 1]
        e0 = errors[change_idx]
        e1 = errors[change_idx + 1]

        denom = e1 - e0
        # 분모가 0인 경우 대비
        mask = denom != np.float64(0.0)

        # 선형보간 근
        x_root_lin = x0 - e0 * (x1 - x0) / denom
        # 분모 0이면 중간값 사용
        x_root_mid = (x0 + x1) / np.float64(2.0)

        x_roots = np.where(mask, x_root_lin, x_root_mid)
        x_cross.extend(x_roots.tolist())

    return x_cross
    # for i in range(len(x_points) - 1):
    #     if errors[i] == 0: # 정확히 0
    #         x_cross.append(x_points[i])
    #     elif signs[i] * signs[i + 1] < 0:
    #         x0, x1 = x_points[i], x_points[i + 1]
    #         err0, err1 = errors[i], errors[i + 1]

    #         if err1 - err0 != 0:
    #             x_root = x0 - err0 * (x1 - x0) / (err1 - err0)
    #         else:
    #             x_root = (x0 + x1) / 2
    #         x_cross.append(x_root)
    # return x_cross

# 구간 변경
def slice_interval(approx_mode: str, intervals: list) -> list:
    if approx_mode == "all":
        return intervals
    
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
        print("singular matrix error!")
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
        # local_x = np.linspace(start, end, 1000)
        local_x = generate_points(start, end)
        local_y = [error_abs(coeff, x, evalF(x)) for x in local_x]
        max_index = np.argmax(local_y)
        max_x = local_x[max_index]
        max_y = local_y[max_index]
        max_point_x.append(max_x)
        max_point_y.append(max_y)
    
    return max_point_x, max_point_y

# 종료조건 판단
def decide_exit(errors: list, threshold: float, print_mode: str) -> bool:
    max_err, min_err = max(errors), min(errors)
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
def calculate_next_remez(coeff, evalF, eb, intervals, err_mode: int):
    ni = []
    max_err = 0

    for start, end in intervals:
        if start == end:
            val = evalF(start) - evalP(coeff, start)
            err = eb.cal_bound(start, coeff, err_mode)/eb.scale
            val1 = val + err
            val2 = val - err
            ni.append([val2, val1])
            max_err = max(val1, max_err)
        else:
            points = generate_points(start, end)
            # length = abs(end - start)
            # # 적절한 간격으로 샘플링, 오차 계산
            # num_points = int(np.clip(np.ceil(length / 1e-3), 50, 2000))
            # points = np.linspace(start, end, num_points)
            
            p_vals1      = np.fromiter((evalP(coeff, p) + eb.cal_bound(p, coeff, err_mode)/eb.scale for p in points), dtype=float, count=points.size)
            p_vals2      = np.fromiter((evalP(coeff, p) - eb.cal_bound(p, coeff, err_mode)/eb.scale for p in points), dtype=float, count=points.size)
            f_vals      = np.fromiter((evalF(p) for p in points), dtype=float, count=points.size)

            err_vals = np.concatenate((np.abs(f_vals - p_vals1), np.abs(f_vals - p_vals2)))
            vals = np.concatenate((p_vals1, p_vals2))
            ni.append([np.min(vals), np.max(vals)])
            
            # 최대오차
            max_err = max(np.max(err_vals), max_err)
    
    # 정렬 - 구간이 긴 쪽이 다음 근사에서 0으로 수렴해야함
    if (ni[0][1]-ni[0][0]) < (ni[1][1]-ni[1][0]):
        ni.reverse()

    # if ni[0][0] > ni[1][0]:
    #     ni.reverse()

    return max_err, ni




