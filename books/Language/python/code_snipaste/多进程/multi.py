from typing import Any, Callable, Iterable, List, Optional
import multiprocessing
import os

T = type('T')

def submit_process_task(task: Callable[..., T],
                        params_list: Iterable[Iterable[Any]],
                        concurrency_level: Optional[int] = None
                        ) -> List[T]:
    pool = multiprocessing.Pool(processes=concurrency_level)
    r = pool.starmap(task, params_list)
    pool.close()
    pool.join()
    return r


def add(x, y):
    return x+y
