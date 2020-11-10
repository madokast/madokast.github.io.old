# a = 1/0 # ZeroDivisionError: division by zero

def div(a,b):
    return a/b

try:
    a = div(1,0)
except ZeroDivisionError as e:
    print(e)
else:
    print("没有异常执行")
finally:
    print("无论有没有异常都执行")



try:
    a = div(1,1)
except ZeroDivisionError as e:
    print(e)
else:
    print("没有异常执行")
finally:
    print("无论有没有异常都执行")

print('-----------------------')


try:
    raise Exception('msg')
except ZeroDivisionError as e:
    print(e)


