"""
307. 区域和检索 - 数组可修改
给定一个整数数组  nums，求出数组从索引 i 到 j  (i ≤ j) 范围内元素的总和，包含 i,  j 两点。

update(i, val) 函数可以通过将下标为 i 的数值更新为 val，从而对数列进行修改。

示例:

Given nums = [1, 3, 5]

sumRange(0, 2) -> 9
update(1, 2)
sumRange(0, 2) -> 8
说明:

数组仅可以在 update 函数下进行修改。
你可以假设 update 函数与 sumRange 函数的调用次数是均匀分布的。
"""

from typing import List


class NumArray:

    def __init__(self, nums: List[int]):
        self.length = len(nums)
        extern_length = 1
        while extern_length < self.length:
            extern_length *= 2

        self.extern_length = extern_length

        self.data_length = self.extern_length*2-1

        self.offset = self.extern_length-1

        self.data = [0]*self.offset + nums + [0]*(self.extern_length-self.length)

        for i in range(self.offset-1, -1, -1):
            self.data[i] = self.data[2*i+1] + self.data[2*i+2]

    def update(self, i: int, val: int) -> None:
        i+=self.offset
        self.data[i] = val
        while i>0:
            i = (i-1)//2
            self.data[i] = self.data[2*i+1] + self.data[2*i+2]
        
        if self.data_length>=3:
            self.data[0] = self.data[1] + self.data[2]

    def sumRange(self, i: int, j: int) -> int:
        sum = 0

        i+=self.offset
        j+=self.offset

        while i<j:
            if i%2==0:
                sum+=self.data[i]
                i+=1
            if j%2==1:
                sum+=self.data[j]
                j-=1
            
            i = (i-1)//2
            j = (j-1)//2
        
        if i==j:
            sum+=self.data[i]
        
        return sum

    def __str__(self) -> str:
        return str(self.data)


# Your NumArray object will be instantiated and called as such:
# obj = NumArray(nums)
# obj.update(i,val)
# param_2 = obj.sumRange(i,j)
if __name__ == "__main__":
    na = NumArray([0,9,5,7,3])
    print(na)

    print(na.sumRange(1,2))

    na.update(1,7)

    print(na.sumRange(1,2))



