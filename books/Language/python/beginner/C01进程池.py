import multiprocessing
import time
import atexit  # 注册程序退出时的钩子函数
from typing import Any, Callable, List, Optional, Type
import functools

T = type('T')


class ComputeAccelerator:
    """
    计算加速器
    """
    __PROCESSES_POOL: multiprocessing.Pool = None

    class ProcessTask:
        """
        描述一个进程任务
        因为 python 具有全局解释器锁，所以 CPU 密集任务无法使用线程加速
        see https://www.cnblogs.com/dragon-123/p/10247252.html
        """

        def __init__(self, *args, **kwargs) -> None:
            """
            args kw 任务参数
            """
            self.args = args
            self.kwargs = kwargs

        def __str__(self) -> str:
            return f"args={self.args}, kwargs={self.kwargs}"

        def __repr__(self) -> str:
            return self.__str__()

    @classmethod
    def submit(cls, func: Callable[['ComputeAccelerator.ProcessTask'], T], tasks: List['ComputeAccelerator.ProcessTask']) -> List[T]:
        if cls.__PROCESSES_POOL is None:
            cls.__init()

        r = cls.__PROCESSES_POOL.map(func, tasks)  # 阻塞

        return r

    @classmethod
    def set_concurrency_level(cls, level: int):
        if cls.__PROCESSES_POOL is not None:
            raise RuntimeError("请在提交任务前设置并发等级")
        if level <= 0:
            raise ValueError("并发等级必须大于0")
        cls.__init(level)

    @classmethod
    def __init(cls, concurrency_level: Optional[int] = None):
        if cls.__PROCESSES_POOL is not None:
            raise RuntimeError("线程池正在运行，不能重复初始化")
        if concurrency_level is None:
            cls.__PROCESSES_POOL: multiprocessing.Pool = multiprocessing.Pool()
        else:
            cls.__PROCESSES_POOL: multiprocessing.Pool = multiprocessing.Pool(
                concurrency_level)

        # 退出函数
        def close():
            print("关闭进程池")
            cls.__PROCESSES_POOL.close()
            cls.__PROCESSES_POOL.terminate()
        # 注册退出函数
        atexit.register(close)


if __name__ == "__main__":
    s = time.time()

    def add4task(task: ComputeAccelerator.ProcessTask):
        print("task.args", task.args)
        add(*task.args, **task.kwargs)

    jobs = [
        ComputeAccelerator.ProcessTask(1, 2, 3),
        ComputeAccelerator.ProcessTask(2, 2, 3),
        ComputeAccelerator.ProcessTask(3, 2, 3),
        ComputeAccelerator.ProcessTask(4, 2, 3),
        ComputeAccelerator.ProcessTask(4, 2, 3, key='a'),
        ComputeAccelerator.ProcessTask(4, 2, 3, key='cc'),
    ]

    print(jobs)

    results = ComputeAccelerator.submit(add4task, jobs)

    print(f"r={results}, time={time.time()-s}")


def add(x, y, z, key='kk'):
    print(f"add,{x},{y},{z},{key}")
    return str(x+y+z)*key