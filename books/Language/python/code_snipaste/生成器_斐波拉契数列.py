def fib(length=5):
    n = 0
    pre = 0
    cur = 1
    while n < length:
        pre, cur = cur, cur + pre
        n += 1
        yield pre


if __name__ == "__main__":
    for i in fib():
        print(i, end=" ")
    print()
    print(list(fib(10)))