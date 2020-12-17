package books.algorithm.java.leetcode.L300;

/**
 * 累加数是一个字符串，组成它的数字可以形成累加序列。
 * 
 * 一个有效的累加序列必须至少包含 3 个数。除了最开始的两个数以外，字符串中的其他数都等于它之前两个数相加的和。
 * 
 * 给定一个只包含数字 '0'-'9' 的字符串，编写一个算法来判断给定输入是否是累加数。
 * 
 * 说明: 累加序列里的数不会以 0 开头，所以不会出现 1, 2, 03 或者 1, 02, 3 的情况。
 * 
 * 示例 1:
 * 
 * 输入: "112358" 输出: true 解释: 累加序列为: 1, 1, 2, 3, 5, 8 。1 + 1 = 2, 1 + 2 = 3, 2 +
 * 3 = 5, 3 + 5 = 8 示例 2:
 * 
 * 输入: "199100199" 输出: true 解释: 累加序列为: 1, 99, 100, 199。1 + 99 = 100, 99 + 100 =
 * 199 进阶: 你如何处理一个溢出的过大的整数输入?
 * 
 * 来源：力扣（LeetCode） 链接：https://leetcode-cn.com/problems/additive-number
 * 著作权归领扣网络所有。商业转载请联系官方授权，非商业转载请注明出处。
 */

public class L305AdditiveNumber {
    public static void main(String[] args) {
        L305AdditiveNumber l = new L305AdditiveNumber();
        System.out.println(l.isAdditiveNumber("000"));
    }

    public boolean isAdditiveNumber(String num) {
        System.out.println(num);
        int len = num.length();
        if (len == 0)
            return true;
        // if (num.charAt(0) == '0')
        // return false;
        boolean ret = false;
        for (int i = 1; i <= (len - 1) / 2; i++) {
            String firstNumberStr = num.substring(0, i);
            long fisrtNumber = Long.parseLong(firstNumberStr);
            if (firstNumberStr.charAt(0) == '0' && firstNumberStr.length() > 1)
                continue;
            for (int j = 1; j <= (len - i) / 2; j++) {
                String secondNumberStr = num.substring(i, i + j);
                if (secondNumberStr.charAt(0) == '0' && secondNumberStr.length() > 1)
                    continue;
                long secondNumber = Long.parseLong(secondNumberStr);
                long sum = fisrtNumber + secondNumber;
                int digitNumber = digitNumber(sum);
                if (i + j + digitNumber > len || i + j >= len)
                    break;
                String thirdNumberStr = num.substring(i + j, i + j + digitNumber);
                long thirdNumber = Long.parseLong(thirdNumberStr);
                if (sum == thirdNumber) {
                    if (i + j + digitNumber == len)
                        return true;
                    else {
                        ret |= isAdditiveNumber(num.substring(i));
                    }
                }
            }
        }

        return ret;
    }

    private int digitNumber(long a) {
        if (a == 0)
            return 1;

        return ((int) Math.log10(a)) + 1;
    }

}
