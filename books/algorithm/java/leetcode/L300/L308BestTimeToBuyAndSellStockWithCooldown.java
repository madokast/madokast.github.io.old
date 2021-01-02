package books.algorithm.java.leetcode.L300;

/**
 * 给定一个整数数组，其中第 i 个元素代表了第 i 天的股票价格 。​
 * 
 * 设计一个算法计算出最大利润。在满足以下约束条件下，你可以尽可能地完成更多的交易（多次买卖一支股票）:
 * 
 * 你不能同时参与多笔交易（你必须在再次购买前出售掉之前的股票）。 卖出股票后，你无法在第二天买入股票 (即冷冻期为 1 天)。 示例:
 * 
 * 输入: [1,2,3,0,2] 输出: 3 解释: 对应的交易状态为: [买入, 卖出, 冷冻期, 买入, 卖出]
 * 
 * 来源：力扣（LeetCode）
 * 链接：https://leetcode-cn.com/problems/best-time-to-buy-and-sell-stock-with-cooldown
 * 著作权归领扣网络所有。商业转载请联系官方授权，非商业转载请注明出处。
 * 
 * 居然一次性做出来了！！！！！！！！！！！！！！ 执行结果： 通过 显示详情 执行用时： 2 ms , 在所有 Java 提交中击败了 17.77% 的用户
 * 内存消耗： 37.7 MB , 在所有 Java 提交中击败了 47.05% 的用户
 * 
 * 提交时间 提交结果 运行时间 内存消耗 语言 
 * 几秒前 通过 2 ms 37.7 MB Java
 */

public class L308BestTimeToBuyAndSellStockWithCooldown {
    public static void main(String[] args) {
        L308BestTimeToBuyAndSellStockWithCooldown l = new L308BestTimeToBuyAndSellStockWithCooldown();
        System.out.println(l.maxProfit(new int[] { 1, 2, 3, 0, 2 }));
    }

    final int STATE_COLD = 0; // 冷冻
    final int STATE_EMPTY = 1; // 空窗 初始
    final int STATE_BUY = 2; // 买入
    final int STATE_HOLD = 3; // 持有
    final int STATE_SELL = 4; // 卖出

    final int LENGTH_OF_STATE = 5;

    final int MIN = Integer.MIN_VALUE;

    public int maxProfit(int[] prices) {
        if (prices == null)
            return 0;
        int len = prices.length;
        if (len <= 1)
            return 0;

        int[][] dp = new int[len][LENGTH_OF_STATE];

        dp[0][STATE_COLD] = MIN; // 不可能
        dp[0][STATE_EMPTY] = 0;
        dp[0][STATE_BUY] = -prices[0];
        dp[0][STATE_HOLD] = MIN;
        dp[0][STATE_SELL] = MIN;

        int ret = 0;

        for (int i = 1; i < len; i++) {
            int p = prices[i];

            dp[i][STATE_COLD] = dp[i - 1][STATE_SELL];
            dp[i][STATE_EMPTY] = Math.max(dp[i - 1][STATE_EMPTY], dp[i - 1][STATE_COLD]);
            dp[i][STATE_BUY] = Math.max(dp[i - 1][STATE_EMPTY], dp[i - 1][STATE_COLD]) - p;
            dp[i][STATE_HOLD] = Math.max(dp[i - 1][STATE_HOLD], dp[i - 1][STATE_BUY]);
            dp[i][STATE_SELL] = Math.max(dp[i - 1][STATE_HOLD], dp[i - 1][STATE_BUY]) + p;

            ret = Math.max(ret, dp[i][STATE_SELL]);
        }

        return ret;
    }
}
