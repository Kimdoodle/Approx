import re
from pathlib import Path
from typing import List, Tuple

# "<계수>x^<지수>" 패턴(공백 허용)
TERM_RE = re.compile(r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*x\^(\-?\d+)', re.I)

def sci_to_desmos(num_str: str) -> str:
    """과학적 표기 'aEb' -> 'a*10^{b}'"""
    s = num_str.strip()
    if 'e' in s.lower():
        base, exp = re.split(r'[eE]', s)
        # base 앞의 +는 제거(가독성), -는 유지하지 않음(부호는 바깥에서 처리)
        neg = base.startswith('-')
        base = base.lstrip('+-')
        exp_int = int(exp)
        mag = f"{base}*10^{{{exp_int}}}"
        # 부호는 바깥에서 처리하므로 크기(magnitude)만 반환
        return mag
    # 부호는 바깥에서 처리하므로 크기만 반환
    return s.lstrip('+-')

def split_sign(s: str) -> Tuple[str, str]:
    """문자열에서 부호와 절댓값 분리. ('+'|'-', magnitude_without_sign)"""
    s = s.strip()
    if s.startswith('-'):
        return '-', s[1:]
    if s.startswith('+'):
        return '+', s[1:]
    return '+', s

def format_magnitude(mag_str: str, power: int) -> str:
    """절댓값 크기와 차수로 항 표현(부호 제외). 지수는 ^{...}."""
    m = sci_to_desmos(mag_str)
    if power == 0:
        return f"{m}"
    elif power == 1:
        return f"{m} x"          # 공백 곱셈
    else:
        return f"{m} x^{{{power}}}"

def parse_polynomial_line(line: str) -> List[Tuple[int, str]]:
    """한 줄에서 (power, coeff_str) 목록 추출."""
    s = line.strip().replace(" ", "")
    matches = TERM_RE.findall(s)
    if not matches:
        raise ValueError(f"다항식 형식을 파싱하지 못했습니다: {line}")
    terms = [(int(p), c) for (c, p) in matches]
    return terms  # [(power, coeff_str)]

def polynomial_to_desmos(terms: List[Tuple[int, str]]) -> str:
    """(power, coeff_str) 리스트를 올바른 부호/연산자로 결합한 다항식 문자열로 변환."""
    if not terms:
        return "0"
    terms_sorted = sorted(terms, key=lambda t: t[0], reverse=True)

    parts = []
    for idx, (p, cstr) in enumerate(terms_sorted):
        sign, mag_raw = split_sign(cstr)
        mag_expr = format_magnitude(mag_raw, p)

        if idx == 0:
            # 첫 항: 음수면 앞에 '-' 추가, 양수면 그대로
            if sign == '-':
                parts.append(f"- {mag_expr}")
            else:
                parts.append(mag_expr)
        else:
            # 이후 항: ' + ' 또는 ' - ' 연산자 명시
            op = '+' if sign == '+' else '-'
            parts.append(f"{op} {mag_expr}")

    return " ".join(parts)

def build_desmos_script(input_path: str, output_path: str = None) -> str:
    """
    입력 파일을 파싱해 Desmos 스크립트 반환:
      f_{i}(x) = <다항식>
    마지막 줄:
      F(x) = f_{n}(...f_{1}(x)...)
    """
    p = Path(input_path)
    lines = [ln for ln in p.read_text(encoding="utf-8").splitlines() if ln.strip()]

    func_lines = []
    for i, raw in enumerate(lines, 1):
        terms = parse_polynomial_line(raw)
        poly_str = polynomial_to_desmos(terms)
        func_lines.append(f"f_{{{i}}}(x) = {poly_str}")

    # 합성 F(x)
    inner = "x"
    for i in range(1, len(func_lines) + 1):
        inner = f"f_{{{i}}}({inner})"
    F_line = f"F(x) = {inner}"

    script = "\n".join(func_lines + [F_line])

    # if output_path:
    #     Path(output_path).write_text(script, encoding="utf-8")

    return script


if __name__ == "__main__":
    p = 2
    e = 15
    in_file = f"doc/coeff_{p}_{e}.txt"           # 입력 파일 경로
    out_file = f"coeff_{p}_{e}_out.txt"     # 출력 파일 경로
    print(build_desmos_script(in_file, out_file))
