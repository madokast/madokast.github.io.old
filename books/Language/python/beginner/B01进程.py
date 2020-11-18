from multiprocessing import Process
import multiprocessing
import time
import sys


def task1():
    print("任务一执行")
    for i in range(1, 10):
        time.sleep(0.1)
        print(f"任务一继续运行{i}")
    


def task2():
    print("任务二执行")
    for i in range(1, 10):
        time.sleep(0.1)
        print(f"任务二继续运行{i}")
    


if __name__ == "__main__":
    p = Process(target=task2)
    p.start()
    task1()


#####################################################


def send(queue: multiprocessing.Queue):
    for i in range(1, 5):
        queue.put(i)
        print(f"send {i} at {i}")
        time.sleep(0.1)


def get(queue: multiprocessing.Queue):
    for i in range(1, 5):
        g = queue.get()
        print(f"get {g} at {i}")


if __name__ == "__main__":
    queue = multiprocessing.Queue(5)
    p1 = multiprocessing.Process(target=send, args=(queue,))
    p2 = multiprocessing.Process(target=get, args=(queue,))

    p1.start()
    p2.start()