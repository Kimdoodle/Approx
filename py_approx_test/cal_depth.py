from math import log2, ceil, floor

class CalData:
    def __init__(self, title="", iter=0):
        self.title = title
        self.iter = iter
        self.depth = 0
        self.add_count = 0
        self.pt_mult_count = 0
        self.ct_mult_count = 0

    def add(self, c2: 'CalData', mode="add"):
        if mode == "add":
            self.depth += c2.depth
        elif mode == "compare":
            self.depth = self.depth if self.depth >= c2.depth else c2.depth
        self.add_count += c2.add_count
        self.pt_mult_count += c2.pt_mult_count
        self.ct_mult_count += c2.ct_mult_count

    def compare(self, other: 'CalData', criteria: str) -> int:
        if criteria == "depth":
            a, b = self.depth, other.depth
        elif criteria == "ctxt":
            a, b = self.ct_mult_count, other.ct_mult_count
        elif criteria == "both":
            ad, bd = self.depth, other.depth
            ac, bc = self.ct_mult_count, other.ct_mult_count
            aa, ba = self.add_count, other.add_count
            
            if ad > bd: return -1
            elif ad < bd: return 1
            else:
                if ac > bc: return -1
                elif ac < bc: return 1
                else:
                    if aa > ba: return -1
                    elif aa < ba: return 1
                    else:
                        return 0

        if a < b:
            return 1
        elif a == b:
            return 0
        elif a > b:
            return -1
        
    def print_params(self, title=False, iter=False):
        if title:
            print(f"Title:\t\t{self.title}")
        if iter:
            print(f"Iter:\t\t{self.iter}")
        print(f"Depth:\t\t{self.depth}")
        print(f"Add:\t\t{self.add_count}")
        print(f"Ptxt Mult:\t{self.pt_mult_count}")
        print(f"Ctxt mult:\t{self.ct_mult_count}")


# 함수 종류 확인
def detect_coeff_type(coeff: list[float]) -> str:
    has_even = any(abs(c) != 0.0 for i, c in enumerate(coeff) if i % 2 == 0)
    has_odd  = any(abs(c) != 0.0 for i, c in enumerate(coeff) if i % 2 == 1)

    if has_even and not has_odd:
        return "even"
    if has_odd and not has_even:
        return "odd"
    return "all"

# 특정 함수의 반복횟수에 대한 연산량, 깊이 계산
def cal_iter(coeff: list[float], iter: int) -> CalData:
    res = cal_coeff(coeff)
    res.depth *= iter
    res.add_count *= iter
    res.ct_mult_count *= iter
    res.pt_mult_count *= iter

    return res


# 특정 함수의 연산량, 깊이 계산
def cal_coeff(coeff: list[float]) -> CalData:
    # 함수가 홀수/짝수/전체인지 확인
    coeff_type = detect_coeff_type(coeff)
    res = CalData()
    max_deg = len(coeff) - 1

    if coeff_type == "all":
        res.ct_mult_count = max_deg - 1
        res.add_count = max_deg
    elif coeff_type == "even":
        res.ct_mult_count = int(max_deg / 2)
        res.add_count = int(max_deg / 2)

    res.pt_mult_count = 1
    res.depth = res.ct_mult_count + 1
    return res
