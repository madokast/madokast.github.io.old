class Solution:
    def isAdditiveNumber(self, num: str) -> bool:
        def dps(index, last1, last2):
            print(last1,last2)
            if index == len(num):
                return True
            sum = last1+last2
            sum_str = str(sum)
            size = len(str(sum_str))
            if index+size <= len(num) and num[index:index+size] == sum_str:
                return dps(index+size, last2, sum)
            else:
                return False

        for i in range(1, len(num)-1):
            if num[0] == '0' and i > 1:
                break
            last1 = int(num[0:i])
            for j in range(1, len(num)):
                if i+j <= len(num)-i and i+j <= len(num)-j:
                    if num[i] == '0' and j > 1:
                        break
                    last2 = int(num[i:i+j])
                    print('-',last1,last2)
                    if dps(i+j, last1, last2):
                        return True

        return False


if __name__ == "__main__":
    print(Solution().isAdditiveNumber('199001200'))
