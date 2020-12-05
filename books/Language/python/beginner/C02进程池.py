import multiprocessing
import time
from typing import Callable, List, Optional
import numpy as np


def submit_process_task(func: Callable, param_list: List[List],concurrency_level: Optional[int] = None):
    pool = multiprocessing.Pool(concurrency_level)
    r = pool.starmap(func, param_list)
    pool.close()
    pool.join()
    return r


def add(x, y=0, key='a'):
    return str(x+y)+"-" + key

class Add:
    @staticmethod
    def fun(x, y=0, key='a'):
        return str(x+y)+"-" + key


if __name__ == "__main__":
    rs = submit_process_task(add, [
        [1], [2], [3], [1, 11], [1, 123, 'bb']
    ],1)

    print(rs)

    print(submit_process_task(Add.fun,[
        [1],
        [np.array([1,2,3])],
        [np.array([1,2,3]),15]
    ]))