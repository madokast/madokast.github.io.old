"""
缓存
304. 二维区域和检索 - 矩阵不可变
给定一个二维矩阵，计算其子矩形范围内元素的总和，该子矩阵的左上角为 (row1, col1) ，右下角为 (row2, col2)。

Range Sum Query 2D
上图子矩阵左上角 (row1, col1) = (2, 1) ，右下角(row2, col2) = (4, 3)，该子矩形内元素的总和为 8。

示例:

给定 matrix = [
  [3, 0, 1, 4, 2],
  [5, 6, 3, 2, 1],
  [1, 2, 0, 1, 5],
  [4, 1, 0, 1, 7],
  [1, 0, 3, 0, 5]
]

sumRegion(2, 1, 4, 3) -> 8
sumRegion(1, 1, 2, 2) -> 11
sumRegion(1, 2, 2, 4) -> 12
说明:

你可以假设矩阵不可变。
会多次调用 sumRegion 方法。
你可以假设 row1 ≤ row2 且 col1 ≤ col2。

执行用时：
84 ms
, 在所有 Python3 提交中击败了
80.50%
的用户
内存消耗：
16.4 MB
, 在所有 Python3 提交中击败了
36.95%
的用户
"""


from typing import List


class NumMatrix:
    def __init__(self, matrix: List[List[int]]):
        self.empty = False
        if matrix is None:
            self.empty = True
            return
        r = len(matrix)
        if r == 0:
            self.empty = True
            return
        c = len(matrix[0])
        if c == 0:
            self.empty = True
            return

        self.dp = [[0] * (c + 1) for i in range(r + 1)]

        for i in range(r):
            for j in range(c):
                self.dp[i + 1][j + 1] = (
                    self.dp[i][j + 1] + self.dp[i + 1][j] - self.dp[i][j] + matrix[i][j]
                )

    def sumRegion(self, row1: int, col1: int, row2: int, col2: int) -> int:
        if self.empty:
            return 0
        else:
            return (
                self.dp[row2 + 1][col2 + 1]
                - self.dp[row1][col2 + 1]
                - self.dp[row2 + 1][col1]
                + self.dp[row1][col1]
            )


# Your NumMatrix object will be instantiated and called as such:
# obj = NumMatrix(matrix)
# param_1 = obj.sumRegion(row1,col1,row2,col2)