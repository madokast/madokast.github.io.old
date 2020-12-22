dp = [[-1]*101 for _ in range(102) ]

def f(m,n):
    """
    最小面值为 m
    组成 n 元
    """
    if dp[m][n] == -1:
        if m==n:
            dp[m][n] = 1
        elif m>n:
            dp[m][n] = 0
        else:
            dp[m][n] = sum([f(i,n-m) for i in range(m,n-m+1)])
    return dp[m][n]

def solve(n):
    """
    组成 n 元
    """
    return sum([f(i,n) for i in range(1,n+1)])

if __name__ == "__main__":
    for i in range(1,100):
        print(i,solve(i))