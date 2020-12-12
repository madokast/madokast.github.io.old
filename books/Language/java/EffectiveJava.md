# 创建和销毁对象

## 1 用静态工厂方法代替构造器

优点：有名字、不必创建新对象、返回子类

缺点：继承上的困难、难以被发现

典型名字：of ValueOf from instance create getType newType type

## 2 当对象实例化需要很对参数，有些必须，有些可选时，使用构造器模式，即 Builder 静态对象

## 3 单例模式下，私有化构造器 和 枚举类

## 4 工具类（不可实例化的类），使用私有空构造器

## 5 优先考虑依赖注入来引入资源

单例类和工具类都不应使用底层资源

## 6 避免创建不必要的对象

new String("abc")

## 7 数组必须显示指明垃圾 ⭐

用数组实现栈，弹栈就是 retuan arr[--size]，但实际上引用未去除，会导致内存泄漏，应该有 arr[--size]=null 才行

## 8 避免使用 finalize clean

利用 finalize 来释放资源是不对的，因为资源宝贵，而 finalize 方法不知何时才会运行

## 9 优先使用 try-with-rource

这是为了代替 try-finally。很多资源需要关闭，如果资源使用中出现异常，很可能进入 finally 中的 close() 方法也会出现异常，这样就覆盖了前面的异常，导致排查困难。


# 对象的通用方法

## 10 覆盖 equals 方法的约定

- 添加了成员 field 的继承和 equals 方法存在本质的冲突？对的

## 11 覆盖 equals 时总要覆盖 hashCode

## 12 始终覆盖 toString

## 13 谨慎使用 clone

- Cloneable 接口的作用是什么？

答：改变父类 clone 方法的行为，如果子类没有实现 Cloneable，super.clone 抛出异常

- 数组实现了 clone 方法？对的

- 比 colne 更好的是使用拷贝构造器

## 14 考虑实现 comparable 接口

- 比较大小时，考虑使用 Type.compare() 更好？对的
