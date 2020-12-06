import multiprocessing
import time
import os
import sys


CONCURRENCY_LOCK_FILE = 'concurrency.lock'


def submit(task, params_list):
    if os.path.exists(CONCURRENCY_LOCK_FILE):
        print(f"{__name__}:已存在锁，不能再创建子进程了")
    else:
        print(f"{__name__}:创建新线程")
        f = open(CONCURRENCY_LOCK_FILE, mode='w')
        f.close()
        p = multiprocessing.Pool()
        r = p.starmap(task, params_list)
        p.close()
        p.join()
        os.remove(CONCURRENCY_LOCK_FILE)
        return r


def add(x, y):
    print(f"name={__name__},x={x},y={y}")
    return x+y

if os.path.exists(CONCURRENCY_LOCK_FILE):
    print(f"name=={__name__},pid={os.getpid()}")
    # sys.exit()

# if __name__ == "__main__":
print(submit(add, [
    [a, b] for a in range(4) for b in range(4)
]))
