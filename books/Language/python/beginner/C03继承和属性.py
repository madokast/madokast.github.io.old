class Father:
    def __init__(self) -> None:
        self.a = 1


class Son(Father):
    def __init__(self) -> None:
        super().__init__()
        self.b = 1111



if __name__ == "__main__":
    f = Father()
    print(f.a)

    s = Son()
    print(s.a)
    print(s.b)