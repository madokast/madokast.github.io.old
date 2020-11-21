from typing import List

"""
301. 删除无效的括号
删除最小数量的无效括号，使得输入的字符串有效，返回所有可能的结果。

说明: 输入可能包含了除 ( 和 ) 以外的字符。

示例 1:

输入: "()())()"
输出: ["()()()", "(())()"]
示例 2:

输入: "(a)())()"
输出: ["(a)()()", "(a())()"]
示例 3:

输入: ")("
输出: [""]
"""

class Solution:
    def removeInvalidParentheses(self, s: str) -> List[str]:
        lefts, rights = 0, 0
        for ch in s:
            if ch == "(":
                lefts += 1
            elif ch == ")":
                if lefts > 0:
                    lefts -= 1
                else:
                    rights += 1

        print(f"要删除{lefts}个左括号，{rights}个右括号")

        ret = set()
        removed = [False] * len(s)
        self.dps(s, removed, 0, ret, lefts, rights)

        return list(ret)

    def dps(self, s, removed, i, ret, lele, leri):
        if self.validate(s, removed):
            ret.add(self.make(s, removed))
        elif i < len(s):
            if s[i] == "(" and lele > 0:
                removed[i] = True
                self.dps(s, removed, i + 1, ret, lele - 1, leri)
                removed[i] = False
            if s[i] == ")" and leri > 0:
                removed[i] = True
                self.dps(s, removed, i + 1, ret, lele, leri - 1)
                removed[i] = False
            self.dps(s, removed, i + 1, ret, lele, leri)

    def make(self, s, removed):
        ret = ""
        for i in range(len(s)):
            if not removed[i]:
                ret += s[i]

        return ret

    def validate(self, s, removed):
        le = 0
        for i in range(len(s)):
            if not removed[i]:
                if s[i] == "(":
                    le += 1
                elif s[i] == ")":
                    if le == 0:
                        return False
                    else:
                        le -= 1

        return le == 0


if __name__ == "__main__":
    s = Solution()
    print(s.removeInvalidParentheses("()())()"))
