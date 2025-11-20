# util_MT.py
import threading
from concurrent.futures import ProcessPoolExecutor
from src_remez.print import write_params

_EXEC = None
_PENDING = 0
_PENDING_LOCK = threading.Lock()
_ALL_DONE = threading.Event()

# ========= 프로세스 풀 제어 =========
_DEFAULT_WORKERS = 32
_MAX_WORKERS: int | None = None  # None이면 기본값(_DEFAULT_WORKERS) 사용


def set_max_workers(num_cores: int | None):
    """
    멀티프로세스에서 사용할 코어(워커) 수를 설정한다.
    num_cores가 None이거나 1 이하이면 기본값(_DEFAULT_WORKERS)을 사용한다.
    이미 실행 중인 풀(_EXEC)이 있다면 종료 후 다시 만든다.
    """
    global _MAX_WORKERS, _EXEC
    if not num_cores or num_cores <= 0:
        _MAX_WORKERS = None
    else:
        _MAX_WORKERS = int(num_cores)

    # 이미 생성된 풀을 재설정
    if _EXEC is not None:
        _EXEC.shutdown(wait=True)
        _EXEC = None


def _get_max_workers() -> int:
    """
    현재 설정된 워커 수를 반환한다.
    별도의 설정이 없으면 기존 코드와의 호환을 위해 32를 사용한다.
    """
    return _MAX_WORKERS if _MAX_WORKERS is not None else _DEFAULT_WORKERS


def _inc():
    global _PENDING
    with _PENDING_LOCK:
        _PENDING += 1
        _ALL_DONE.clear()


def _dec():
    global _PENDING
    with _PENDING_LOCK:
        _PENDING -= 1
        if _PENDING == 0:
            _ALL_DONE.set()


def _ensure_exec(max_workers: int | None = None):
    """
    멀티프로세스 풀 초기화 보장.
    max_workers를 명시하면 해당 값으로, 생략하면 set_max_workers로 설정된 값
    또는 기본값(_DEFAULT_WORKERS)으로 풀을 생성한다.
    """
    global _EXEC
    if _EXEC is None:
        if max_workers is None:
            max_workers = _get_max_workers()
        _EXEC = ProcessPoolExecutor(max_workers=max_workers)


def wait_for_all():
    """
    모든 제출 작업 완료 대기
    """
    _ALL_DONE.wait()


def shutdown_exec():
    """
    프로세스 풀 종료
    """
    global _EXEC
    if _EXEC is not None:
        _EXEC.shutdown(wait=True)
        _EXEC = None

import multiprocessing
from multiprocessing.managers import SyncManager

# Manager 및 공유 객체들
_MANAGER: SyncManager | None = None
_STORE = None          # { depth: remezData }
_PRE_LOCKS = None    # { pre: RLock }
_GLOBAL_LOCK = None    # precision용 락 생성 시에만 잠그는 용도

def shared_clear():
    # 공유 메모리 초기화
    _ensure_manager()
    with _GLOBAL_LOCK:
        _STORE.clear()
        _PRE_LOCKS.clear()
        
def _ensure_manager():
    global _MANAGER, _STORE, _PRE_LOCKS, _GLOBAL_LOCK
    if _MANAGER is None:
        _MANAGER = multiprocessing.Manager()
        _STORE = _MANAGER.dict()
        _PRE_LOCKS = _MANAGER.dict()
        _GLOBAL_LOCK = _MANAGER.RLock()

def _get_pre_lock(pre: int):
    """
    pre별 전용 락을 반환. 없으면 생성.
    """
    _ensure_manager()
    with _GLOBAL_LOCK:   # 락 생성/등록 시에만 전역 잠금
        lock = _PRE_LOCKS.get(pre)
        if lock is None:
            lock = _MANAGER.RLock()
            _PRE_LOCKS[pre] = lock
        return lock
    
def shared_init():
    """
    공유 스토어 초기화 (멀티프로세스 전역)
    """
    _ensure_manager()

def shared_get(key: int):
    _ensure_manager()
    with _STORE_LOCK:
        return _STORE.get(key, None)

        
def shared_update_if_better(p_num: int, newData):
    """
    compare == True인 경우에만 교체.
    기존 파일 로직(기존.compare(new) == False → 저장 안 함)을
    공유 메모리로 치환.
    """
    _ensure_manager()
    with _STORE_LOCK:
        prev = _STORE.get(p_num, None)
        if prev is None:
            _STORE[p_num] = newData
            return True

        # prev.compare(newData) == True 이면 newData가 더 좋다고 판단하고 교체
        try:
            res = prev.compare(newData)
        except Exception:
            # 비교 실패 시 안전하게 교체하지 않음
            return False

        if res is True:
            _STORE[p_num] = newData
            write_params(p_num, 30, newData)
            print(f"depth = {newData.total_CalData.depth}")
            return True
        else:
            return False