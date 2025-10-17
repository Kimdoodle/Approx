from .basic import evalP, error_func, error_abs
import numpy as np
from .print import debug_print

# u(x)연산
def evalU(coeff, x, evalF, h=1e-10) -> int:
    xp = x + h
    xm = x - h

    # 경계값 고려
    if evalF(xp) is False:
        xp = x - h
    if evalF(xm) is False:
        xm = x + h  

    try:
        temp1 = evalP(coeff, xp) - evalF(xp)
        temp2 = evalP(coeff, x) - evalF(x)
        temp3 = evalP(coeff, xm) - evalF(xm)
    except TypeError as e:
        print(f"xp: {xp}, x: {x}, xm: {xm}")
        print(f"P(x): {evalP(coeff, xp)}, {evalP(coeff, x)}, {evalP(coeff, xm)}")
        print(f"f(x): {evalF(xp)}, {evalF(x)}, {evalF(xm)}")
        
    dd = (temp1 - 2*temp2 + temp3) / (h*h)
    # print(f"x: {x}\ntemp1:{temeps1}, temp2:{temp2}, temp3:{temp3}\nmu: {dd}")
    eps = 1e-5
    if   dd > eps: return  1
    elif dd < -eps: return -1
    else:           return  0
    
# 극점 조건 확인
def check_multi_interval_conditions(points: list, E: float, coeff, evalF) \
    -> tuple[float, bool]:
    
    # 1번 조건 - Local Extreme Value Conditon
    cond1_list = [evalU(coeff, p, evalF) * error_func(coeff, p, evalF(p)) for p in points]
    if min(cond1_list) < E:
        return 0, False
    
    # 2번 조건 - Alternating condition
    prev = evalU(coeff, points[0], evalF)
    for i in range(1, len(points)):
        next = evalU(coeff, points[i], evalF)
        if prev * next != -1:
            return 0, False
        prev = next

    # 3번 조건 - Maximum absolute sum condition
    cond3_list = [error_abs(coeff, p, evalF(p)) for p in points]
    sum_value = sum(cond3_list)
    return sum_value, True


def find_best_combination(points_x: list, points_y: list, number: int, coeff, evalF, print_mode) -> tuple:
    debug_print("Finding best combination...", print_mode)
    orig_x, orig_y = points_x, points_y
    h = 1e-10
    while h <= 1e-1:
        # 복사본 생성
        px = orig_x.copy()
        py = orig_y.copy()
        
        # 1) mu 계산
        mu = [evalU(coeff, x, evalF, h) for x in px]
        # print(f"mu: {mu}")

        # 제거 헬퍼
        def remove(idx):
            px.pop(idx); py.pop(idx); mu.pop(idx)
        
        # 2) 인접한 u(x) 같은 부호 제거
        i = 0
        while i < len(px) - 1:
            if mu[i] * mu[i+1] == 1:
                # 오차(points_y)가 작은 쪽 제거
                rem = i if py[i] < py[i+1] else i+1
                remove(rem)
            else:
                i += 1
        
        # 3) 합 배열 T 생성
        def build_T():
            T = []
            for j in range(len(py) - 1):
                T.append((py[j] + py[j+1], j))
            return sorted(T, key=lambda x: x[0])
        
        # 4) 원하는 개수가 될 때까지 제거
        while len(px) > number:
            T = build_T()
            n = len(px)
            
            # A. 한 개만 넘을 때
            if n == number + 1:
                # 양 끝에서 오차 작은 쪽 제거
                rem = 0 if py[0] < py[-1] else n-1
                remove(rem)
            
            # B. 두 개 넘을 때
            elif n == number + 2:
                end_sum = py[0] + py[-1]
                T.append((end_sum, 'end'))
                T.sort(key=lambda x: x[0])
                
                smallest, idx = T[0]
                if idx == 'end':
                    # 끝 양쪽 모두 제거
                    remove(-1)
                    remove(0)
                else:
                    remove(idx+1)
                    remove(idx)
            
            # C. 그 외
            else:
                smallest, idx = T[0]
                if idx == 0:
                    remove(0)
                elif idx == n-2:
                    remove(-1)
                else:
                    remove(idx+1)
                    remove(idx)
        
        # 정확히 number 개가 되었으면 반환
        if len(px) == number:
            # print(f"Converged with eps = {eps}")
            return px, py
        
        # 아니면 eps 10배 증가시키고 재시도
        debug_print("Decreasing h...", print_mode)
        # print(f"eps={eps} 로는 {len(px)}개, 목표는 {number}개 → eps를 10배로 키워 재시도")
        h *= 10
            
        
        
def find_best_combination_v2(points_x: list, points_y: list, coeff, evalF) -> tuple:
    # 1) sgn 계산
    sgn = []
    for x in points_x:
        diff = evalF(x) - evalP(coeff, x)
        sgn.append(1 if diff > 0 else -1)

    # 2) 연속 구간(run)별로 points_y 최대값 인덱스만 남기기
    #    - 원본 인덱스를 기준으로 보존/제거를 결정
    keep_indices = []

    i = 0
    n = len(points_x)
    while i < n:
        # 현재 run의 시작
        run_sign = sgn[i]
        run_start = i
        # 같은 부호가 이어지는 끝까지 확장
        i += 1
        while i < n and sgn[i] == run_sign:
            i += 1
        run_end = i  # [run_start, run_end) 구간이 동일 부호

        # 이 run에서 points_y가 최대인 원소의 '전역 인덱스'를 선택
        best_idx = run_start
        best_y = points_y[run_start]
        for j in range(run_start + 1, run_end):
            if points_y[j] > best_y:
                best_y = points_y[j]
                best_idx = j

        keep_indices.append(best_idx)

    # 3) 선택된 인덱스만 남겨 필터링(원소 순서는 원본 순서 유지)
    #    keep_indices는 run별 하나씩이므로 이미 오름차순
    filtered_x = [points_x[i] for i in keep_indices]
    filtered_y = [points_y[i] for i in keep_indices]

    return filtered_x, filtered_y


def find_best_combination_v3(points_x: list, points_y: list, extra_number: int, coeff, evalF) -> tuple:
    """
    Remove `extra_number` points from (points_x, points_y) following the rules:
      - Compute error signs err_sgn = sign(evalP(coeff, x) - evalF(x)).
      - If the error signs alternate perfectly:
          * If only one removal is needed and a zero exists, remove the zero.
          * If still need to remove, remove one from the front.
          * If still need to remove, remove one from the back.
          * If more removals are needed, keep removing zeros first (if any), then front/back.
      - Otherwise (not perfectly alternating), remove points with the smallest absolute error first.
    Returns:
      (new_points_x, new_points_y)
    """
    def sgn(x):
        if x > 0: return 1
        elif x < 0: return -1
        else: return 0

    def is_perfect_alternating(signs):
        """Check if non-zero signs alternate perfectly (+1, -1, +1, -1, ...)."""
        nz = [s for s in signs if s != 0]
        if len(nz) < 2:
            return False
        return all(nz[i] * nz[i - 1] == -1 for i in range(1, len(nz)))

    # Local copies so we don't mutate the caller's lists
    xs = list(points_x)
    ys = list(points_y)

    if extra_number <= 0 or len(xs) == 0:
        return xs, ys

    # Initial error signs (fixing the missing parenthesis from the prompt)
    err_sgn = [sgn(evalP(coeff, x) - evalF(x)) for x in xs]

    # Case 1: perfect alternation
    if is_perfect_alternating(err_sgn):
        # If exactly one removal and there exists a zero, remove that zero first
        while extra_number > 0:
            if extra_number == 1 and 0 in err_sgn:
                idx = err_sgn.index(0)
                del err_sgn[idx]
                del xs[idx]
                del ys[idx]
                extra_number -= 1
                break  # done

            # Remove zeros first (one by one) if any remain
            if 0 in err_sgn and extra_number > 0:
                idx = err_sgn.index(0)
                del err_sgn[idx]
                del xs[idx]
                del ys[idx]
                extra_number -= 1
                continue

            # Then remove from the front
            if extra_number > 0 and xs:
                del err_sgn[0]
                del xs[0]
                del ys[0]
                extra_number -= 1

            # Then remove from the back
            if extra_number > 0 and xs:
                del err_sgn[-1]
                del xs[-1]
                del ys[-1]
                extra_number -= 1

        return xs, ys

    # Case 2: not perfectly alternating → remove smallest |error|
    # Compute absolute errors
    errs = [abs(evalP(coeff, x) - evalF(x)) for x in xs]
    # Indices sorted by ascending error
    idx_sorted = sorted(range(len(errs)), key=lambda i: errs[i])
    to_remove = set(idx_sorted[:max(0, min(extra_number, len(xs)))])

    new_xs = [x for i, x in enumerate(xs) if i not in to_remove]
    new_ys = [y for i, y in enumerate(ys) if i not in to_remove]

    return new_xs, new_ys

    

    