def f1():
    f2()


def f2():
    print_traceback()


def print_traceback():
    import sys

    f = sys._getframe()
    while f != None:
        print(f)
        f = f.f_back


if __name__ == "__main__":
    f1()