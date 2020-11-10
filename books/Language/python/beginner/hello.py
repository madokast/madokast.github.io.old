def solve(len,num):
    if num==1:
        return len
    else:
        m = 0
        for i in range(1,len):
            m = max(m,solve(len-i,num-1)*i)
        return m

ans = solve(8,3)
print(ans)

