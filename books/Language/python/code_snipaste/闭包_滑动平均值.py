def make_averager():
    count = 0
    total = 0
    def averager(new_value):
        nonlocal count, total
        count += 1
        total += new_value
        return total / count
    return averager

a = make_averager()
print(a(1))
print(a(2))
print(a(3))
print(a(4))

