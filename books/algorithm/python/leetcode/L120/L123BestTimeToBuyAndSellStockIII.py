from typing import List


class Solution:
    def maxProfit(self, prices: List[int]) -> int:
        if prices is None or len(prices)==0:
            return 0
        b1 = -10000000
        b2 = -10000000
        s1 = 0
        s2 = 0
        for p in prices:
            b1 = max(b1,-p)
            s1 = max(s1,b1+p)
            b2 = max(b2,s1-p)
            s2 = max(s2,b2+p)
        
        return max(s1,s2)

if __name__ == "__main__":
    s = Solution()
    print(s.maxProfit([3,3,5,0,0,3,1,4]))