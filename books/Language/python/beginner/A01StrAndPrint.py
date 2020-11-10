print("打印字面量")

name = "madoka"

print("打印变量" + name)

print("打印多个字符串", "str1", "str2", "默认空格分割")

print("打印多个字符串", "str1", "str2", "默认空格分割")

print("打印多个字符串", "str1", "str2", 'sep="修改分割"', sep="#")

print(r"\n")

# 注释

# 占位符

print("age=%d,name=%s,pi=%5.2f" % (12, "mdk", 3.14))

str1 = "age=%d,name=%s" % (12, "mdk")

print(str1)

# print(1 + "1")

print("age={},name={},pi={}".format(12, "mdk", 3.14))
