import os
import fcntl

from parse import parse_into_remezData, get_file
from src_remez.print import coeff2txt

# ========= 파일 잠금 =========
_LOCKFILE_PATH = None

def init_lockfile(doc_dir: str = "doc"):
    global _LOCKFILE_PATH
    os.makedirs(doc_dir, exist_ok=True)
    _LOCKFILE_PATH = os.path.join(doc_dir, ".coeff_lock")

def _with_file_lock(fn, *args, **kwargs):
    if _LOCKFILE_PATH is None:
        init_lockfile()
    with open(_LOCKFILE_PATH, "w") as lf:
        fcntl.flock(lf, fcntl.LOCK_EX)
        try:
            return fn(*args, **kwargs)
        finally:
            fcntl.flock(lf, fcntl.LOCK_UN)

# ========= 파일 읽기/쓰기 =========
def read_params_threadsafe(p_num):
    return _with_file_lock(get_file, p_num)

def write_params_threadsafe(p_num, newData, doc_dir: str = "doc"):
    return _with_file_lock(write_params, p_num, newData, doc_dir)

def write_params(p_num, e_num, newData, doc_dir: str = "doc"):
    """
    기존 결과와 비교하여 더 나은 결과만 저장.
    """
    # try:
    #     depth = newData.total_CalData.depth
    #     existing_files = get_file(p_num)

    #     if len(existing_files) == 0:
    #         pass
    #     else:
    #         for file in existing_files:
    #             prev_data = parse_into_remezData(file)
    #             if prev_data is not False:
    #                 res = prev_data.compare(newData)
    #                 if res is False:
    #                     # 기존 파일이 더 좋으면 저장하지 않음
    #                     # print(f"False, {prev_data.total_CalData.depth} vs {depth}")
    #                     return
    #             else:
    #                 # 파싱 실패 파일은 제거
    #                 os.remove(file)
    # except Exception as e:
    #     print(e)
    #     return

    # 새 파일 생성
    os.makedirs(doc_dir, exist_ok=True)
    filepath = os.path.join(doc_dir, f"coeff_{p_num}_{e_num}.txt")
    with open(filepath, "w", encoding="utf-8") as f:
        for coeff in newData.coeff_log:
            f.write(coeff2txt(coeff))
            f.write('\n')
