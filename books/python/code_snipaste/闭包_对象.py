def student(name: str, age: int) -> dict:
    def hello() -> None:
        print(f"我是{name}，年龄{age}")

    def setName(newName: str) -> None:
        nonlocal name
        name = newName

    def getName() -> str:
        return name

    return {"hello": hello, "setName": setName, "getName": getName}


s1 = student("mdk", 12)
s2 = student("miao", 13)
print(s1)
print(s2)

s1["hello"]()
s2["hello"]()

print(s1["getName"]())
s1["setName"]("zrx")
s1["hello"]()
print(s1["getName"]())
