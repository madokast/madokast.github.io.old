a = 10

while True:
    try:
        b = a * 5
        if b > 0:
            c = 1 // 0
        else:
            break
    except Exception as e:
        print("异常", e)
        a -= 1


print("继续运行")