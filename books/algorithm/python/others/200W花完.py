if __name__ == "__main__":
    # 200w
    m = 200*100*100
    # 每日花200元
    cost_per_day = 200
    # day
    day = 0
    while m >0:
        interest = m * (0.6/10000)
        m+=interest
        m-=cost_per_day
        day+=1
        print(f"第{day}天，利率{interest}，使用{cost_per_day}后剩余{m}")