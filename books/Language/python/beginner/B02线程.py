import time
import threading


def task1():
    for i in range(5):
        time.sleep(0.1)
        print(f"task1-{i}")


def task2():
    for i in range(5):
        time.sleep(0.1)
        print(f"task2-{i}")


if __name__ == "__main__":
    t1 = threading.Thread(target=task1)
    t2 = threading.Thread(target=task2)

    t1.start()
    t2.start()

    t1.join()
    t2.join()

###################################

global_var = 0


def adder():
    global global_var
    for i in range(1000000):
        global_var+=1
    print(global_var)

if __name__ == "__main__":
    t1 = threading.Thread(target=adder)
    t2 = threading.Thread(target=adder)

    t1.start()
    t2.start()

    t1.join()
    t2.join()
    print(global_var)
