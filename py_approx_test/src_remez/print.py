
def debug_print(message: str, mode: str, end='\n') -> None:
    if mode == "debug":
        print(message, end=end)
        
def print_coeff(coeff: list, pre: int):
    print(f"Final polynomial\t: ", end='')
    for i in range(len(coeff)):
        if abs(coeff[i]) < 10 ** (-pre):
            continue
        if coeff[i] >= 0:
            print(f"+{coeff[i]:.{pre}f}(x^{i})", end="\t")
        else:
            print(f"{coeff[i]:.{pre}f}(x^{i})", end="\t")
    print("")


def print_intervals(intervals: list, pre: int):
    for interval in intervals:
        print(f"[{float(interval[0]):.{pre}f}, {float(interval[1]):.{pre}f}]", end=' ')
    print()

def interval2txt(intervals: list) -> str:
    txt = "["
    txt += f"[{intervals[0][0]}, {intervals[0][1]}]"
    txt += ", "
    txt += f"[{intervals[1][0]}, {intervals[1][1]}]"
    txt += "]"
    return txt

        
def coeff2txt(coeff: list) -> str:
    txt = ""
    for i, ele in enumerate(coeff):
        if ele != 0.0:
            if ele > 0:
                txt += f"+{ele}x^{i}"
            else:
                txt += f"{ele}x^{i}"
    return txt