from math import log2, ceil, floor

class CalData:
    def __init__(self, title="", iter=0):
        self.title = title
        self.iter = iter
        self.depth = 0
        self.cadd = 0
        self.pmult = 0
        self.cmult = 0

    def add(self, c2: 'CalData', mode="add"):
        if mode == "add":
            self.depth += c2.depth
        elif mode == "compare":
            self.depth = self.depth if self.depth >= c2.depth else c2.depth
        self.cadd += c2.cadd
        self.pmult += c2.pmult
        self.cmult += c2.cmult

    def compare(self, other: 'CalData', criteria: str) -> int:
        if criteria == "depth":
            a, b = self.depth, other.depth
        elif criteria == "ctxt":
            a, b = self.cmult, other.cmult
        elif criteria == "both":
            ad, bd = self.depth, other.depth
            ac, bc = self.cmult, other.cmult
            aa, ba = self.cadd, other.cadd
            
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
        print(f"Add:\t\t{self.cadd}")
        print(f"Ptxt Mult:\t{self.pmult}")
        print(f"Ctxt mult:\t{self.cmult}")


# 함수 종류 확인
def check_coeff_type(coeff: list):
    degrees = [i for i, c in enumerate(coeff) if c != 0]
    if not degrees:
        return "etc"

    max_deg = max(degrees)
    
    # 1. 홀수차 함수
    if all(d % 2 == 1 for d in degrees):
        return "odd"

    # 2. 짝수차 함수
    if all(d % 2 == 0 for d in degrees):
        return "even"

    # 3. Cleanse 함수
    if max_deg % 2 == 1 and all(d in {max_deg, max_deg - 1} for d in degrees):
        return "cl"

    # 4. 기타
    return "etc"

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
    res.cadd *= iter
    res.cmult *= iter
    res.pmult *= iter
    return res

# 특정 함수의 연산량, 깊이 계산
def cal_coeff(coeff: list[float]) -> CalData:        
    # 함수가 홀수/짝수/전체인지 확인
    coeff_type = check_coeff_type(coeff)
    res = CalData()
    max_deg = len(coeff) - 1
    data = {
        2: (2, 1, 1, 2),
        3: (2, 2, 2, 3),
        4: (3, 2, 3, 4),
        5: (3, 3, 3, 5),
        6: (3, 4, 3, 6)
    }
    match coeff_type:
        case "cl":
            res.cmult = int((max_deg+1) / 2)
            res.cadd = 1
            res.pmult = 1
            res.depth = int(ceil(log2(max_deg + 1)))
            return res
        case "etc":
            res.depth, res.cmult, res.pmult, res.cadd = data[max_deg]
            # print(res.depth, res.cmult, res.pmult, res.cadd)
            
            # res.ct_mult_count = max_deg - 1
            # res.add_count = max_deg
            # res.depth = int(ceil(log2(max_deg + 1)))
        case "even":
            if max_deg <= 4:
                res.cmult = int(max_deg / 2)
                res.cadd = int(max_deg / 2)
                res.depth = int(ceil(log2(max_deg + 1)))
            else:
                res.depth = 3
                res.cmult = 3
                res.pmult = 2
                res.cadd = 3
                return res
                
        case "odd":
            res.cmult = int(max_deg / 2) + 1
            res.cadd = int(max_deg / 2)
            res.depth = res.cmult

    res.pmult = 1

    return res
