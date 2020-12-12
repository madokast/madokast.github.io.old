package books.algorithm.java.leetcode.L300;

/**
 * 缓存问题 给定一个整数数组  nums，求出数组从索引 i 到 j（i ≤ j）范围内元素的总和，包含 i、j 两点。
 * 
 * 实现 NumArray 类：
 * 
 * NumArray(int[] nums) 使用数组 nums 初始化对象 int sumRange(int i, int j) 返回数组 nums
 * 从索引 i 到 j（i ≤ j）范围内元素的总和，包含 i、j 两点（也就是 sum(nums[i], nums[i + 1], ... ,
 * nums[j])）  
 * 
 * 示例：
 * 
 * 输入： ["NumArray", "sumRange", "sumRange", "sumRange"] [[[-2, 0, 3, -5, 2,
 * -1]], [0, 2], [2, 5], [0, 5]] 输出： [null, 1, -1, -3]
 * 
 * 解释： NumArray numArray = new NumArray([-2, 0, 3, -5, 2, -1]);
 * numArray.sumRange(0, 2); // return 1 ((-2) + 0 + 3) numArray.sumRange(2, 5);
 * // return -1 (3 + (-5) + 2 + (-1)) numArray.sumRange(0, 5); // return -3
 * ((-2) + 0 + 3 + (-5) + 2 + (-1))
 * 
 * 来源：力扣（LeetCode） 链接：https://leetcode-cn.com/problems/range-sum-query-immutable
 * 著作权归领扣网络所有。商业转载请联系官方授权，非商业转载请注明出处。
 */

public class L0303RangeSumQueryImmutable {
    public static void main(String[] args) {

    }

}

class NumArray {

    int[] data;
    int[] sum2i;

    public NumArray(int[] nums) {
        int len = nums.length;
        data = nums;
        sum2i = new int[len];
        if (len > 0) {
            sum2i[0] = nums[0];
            for (int i = 1; i < len; i++) {
                sum2i[i] = nums[i] + sum2i[i - 1];
            }
        }

    }

    public int sumRange(int i, int j) {
        return sum2i[j] - sum2i[i] + data[i];
    }
}

/**
 * Your NumArray object will be instantiated and called as such: NumArray obj =
 * new NumArray(nums); int param_1 = obj.sumRange(i,j);
 */