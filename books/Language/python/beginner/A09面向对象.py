from typing import overload


print(isinstance(1, object))  # True


class Student:
    pass


print(Student)  # <class '__main__.Student'>

print(Student())  # <__main__.Student object at 0x0000026355591FD0>

print("类的属性")


class Person:
    MAX_AGE = 100  # 类对象
    __MIN_AGE = 0  # 私有方法 type object 'Person' has no attribute '__MIN_AGE'

    def __init__(self, age):  # 构造器。魔术方法之一，魔术方法就是__xx__()
        self.age = min(age, Person.MAX_AGE)  # 实例对象

    def normalMethod(self):  # 普通方法 实例方法
        print(f"普通方法就是实例方法 age={self.age}")

    @classmethod  # 类方法。至少有一个入参，入参为类自身
    def clsMethod(cls, val):
        print(f"类方法{cls.MAX_AGE} val={val}")

    @staticmethod  # 静态方法。可以没有入参
    def staMethod():
        print(f"静态方法{Person.MAX_AGE} {Person.__MIN_AGE}")  # 静态方法可以访问Person的私有对象


print(Person.MAX_AGE)
print(Person(10).MAX_AGE)
print(Person(10).age)
print(Person(11).normalMethod())
print(Person.clsMethod(Person(1)))
print(Person(1000).age)
# print(Person.__MIN_AGE) # type object 'Person' has no attribute '__MIN_AGE'
print(Person.staMethod())

print("--------------")


class Cct:
    def magnetAt(self, p):
        return p + 1


def fun(cct):
    return cct.magnetAt(5)


Cct.magnetAlone = fun

c = Cct()

print(c)

print(c.magnetAt(4))

print(c.magnetAlone())


print("-----------私有化--------------")


class A:
    def __init__(self, name):
        self.__name = name

    def getName(self):
        return self.__name

    def setName(self, name):
        self.__name = name


a = A("zrx")
# print(a.__name)  # 'A' object has no attribute '__name'
print(a.getName())
a.setName("mdk")
print(a.getName())

print(dir(a))
print(a._A__name)


class A:
    def __init__(self, name):
        self.__name = name

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, name):
        self.__name = name


a = A("1111")
print(a.name)
a.name = "22"
print(a.name)


print("-------------继承---------------")


class Father:
    def __init__(self, name: str) -> None:
        self.name = name

    def hello(self) -> None:
        print(f"hello {self.name}")


class Son(Father):
    def __init__(self, name: str, age: int) -> None:
        super(Son, self).__init__(name)  # 相当于Java super(val)
        self.age = age

    def hello(self) -> None:  # 重写方法
        print(f"hello {self.name} {self.age}")

    def fun(self) -> None:
        print("子类的方法")


s = Son("abc", 12)
s.hello()
s.fun()
f = Father("zzz")
f.hello()
print(Son.__mro__)