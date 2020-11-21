package books.algorithm.java.leetcode.L300;

import java.util.Arrays;

/**
 * 给定一个无序的整数数组，找到其中最长上升子序列的长度。
 * 
 * 示例:
 * 
 * 输入: [10,9,2,5,3,7,101,18] 输出: 4 解释: 最长的上升子序列是 [2,3,7,101]，它的长度是 4。 说明:
 * 
 * 可能会有多种最长上升子序列的组合，你只需要输出对应的长度即可。 你算法的时间复杂度应该为 O(n2) 。 进阶:
 * 你能将算法的时间复杂度降低到 O(nlog n) 吗?
 * 
 * 来源：力扣（LeetCode）
 * 链接：https://leetcode-cn.com/problems/longest-increasing-subsequence
 * 著作权归领扣网络所有。商业转载请联系官方授权，非商业转载请注明出处。
 */

public class L0300LongestIncreasingSubsequence {
    public static void main(String[] args) {
        L0300LongestIncreasingSubsequence a = new L0300LongestIncreasingSubsequence();
        System.out.println(a.lengthOfLIS(new int[] { 10, 9, 2, 5, 3, 7, 101, 18 })); // 4
        System.out.println(a.lengthOfLIS(new int[] { 1, 3, 6, 7, 9, 4, 10, 5, 6 })); // 6
    }

    private int lengthOfLIS(int[] nums) {
        if (nums == null || nums.length == 0)
            return 0;
        if (nums.length == 1)
            return 1;

        int len = nums.length;

        int[] dp = new int[len];

        Arrays.fill(dp, 1);

        for (int i = 1; i < len; i++) {
            // dp[i] = ?
            for (int j = i - 1; j >= 0; j--) {
                if (nums[i] > nums[j] && j + 2 > dp[i]) {
                    dp[i] = Math.max(dp[i], dp[j] + 1);
                }
            }
        }

        return max(dp);
    }

    private int max(int[] dp) {
        int m = 0;
        for (int i : dp) {
            m = Math.max(m, i);
        }
        return m;
    }
}
