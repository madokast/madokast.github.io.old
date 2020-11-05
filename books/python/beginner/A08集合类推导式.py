names = ["abc", "qqqqqq", "aa"]

shortNames = [n for n in names if len(n) <= 3]

print(shortNames)

lens = [len(n) for n in names]

print(lens)

print([i for i in range(101) if i % 3 == 0])

print("---------------")

print([(x, y) for x in range(5) if x % 2 == 0 for y in range(5) if y % 2 != 0])

arr1 = ["a", "b", "c"]
arr2 = [1, 2, 3]
cross = [(x, y) for x in arr1 for y in arr2]
print(cross)


m = {10: "aaa", 2: "b"}
print({k: len(v) for k, v in m.items()})
