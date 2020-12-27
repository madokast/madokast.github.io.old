package books.algorithm.java.leetcode.L300;

import java.util.Arrays;

/**
 * 给定一个整数数组  nums，求出数组从索引 i 到 j  (i ≤ j) 范围内元素的总和，包含 i,  j 两点。
 * 
 * update(i, val) 函数可以通过将下标为 i 的数值更新为 val，从而对数列进行修改。
 * 
 * 示例:
 * 
 * Given nums = [1, 3, 5]
 * 
 * sumRange(0, 2) -> 9 update(1, 2) sumRange(0, 2) -> 8 说明:
 * 
 * 数组仅可以在 update 函数下进行修改。 你可以假设 update 函数与 sumRange 函数的调用次数是均匀分布的。
 * 
 * 来源：力扣（LeetCode） 链接：https://leetcode-cn.com/problems/range-sum-query-mutable
 * 著作权归领扣网络所有。商业转载请联系官方授权，非商业转载请注明出处。
 */

/**
 * 这个题目主要学一个叫做线段树的东西
 */

public class L307RangeSumQueryMutable {
    public static void main(String[] args) {
        System.out.println("L307RangeSumQueryMutable");
        int[] nums = {1,2,3,4,5,6};

        L307RangeSumQueryMutable l = new L307RangeSumQueryMutable();

        NumArray na = l.new NumArray(nums);

        System.out.println(na);

        System.out.println(na.sumRange(1, 4));

        na.update(3, 14);

        System.out.println(na.sumRange(1, 4));

        System.out.println("----------------");

        na = l.new NumArray(new int[]{0,9,5,7,3});

//["NumArray","  sumRange","sumRange","sumRange","update","update","update","sumRange","update","sumRange","update"]
//[[[0,9,5,7,3]],[4,4],     [2,4],    [3,3],      [4,5],   [1,7],   [0,8],   [1,2],    [1,9],   [4,4],     [3,4]]

        System.out.println(na);

        na.update(4, 5);
        na.update(1, 7);
        na.update(0, 8);

        System.out.println(na);

        System.out.println(na.sumRange(1, 2));
    }

    class NumArray {

        int[] data;

        int dataLength;

        int offset;

        public NumArray(int[] nums) {
            int len = nums.length;
            int externLength = 1;
            while (externLength < len)
                externLength *= 2;

            int dataLength = externLength * 2 - 1;

            data = new int[dataLength];
            this.dataLength = dataLength;
            this.offset = externLength - 1;

            System.arraycopy(nums, 0, data, externLength - 1, len);

            for (int i = externLength - 2; i >= 0; i--) {
                data[i] = data[2 * i + 1] + data[2 * i + 2];
            }
        }

        public void update(int i, int val) {
            int index = i + offset;
            data[index] = val;

            while (index > 0) {
                index = (index - 1) / 2;
                data[index] = data[index * 2 + 1] + data[index * 2 + 2];
            }

            if(dataLength>=3) data[0] = data[1] + data[2];
        }

        public int sumRange(int i, int j) {
            int sum = 0;

            i += offset;
            j += offset;

            while (i < j) {
                if (i % 2 == 0) {
                    sum += data[i];
                    i++;
                }

                if (j % 2 == 1) {
                    sum += data[j];
                    j--;
                }

                i = (i - 1) / 2;
                j = (j - 1) / 2;
            }

            if(i==j) sum+=data[i];

            return sum;
        }

        @Override
        public String toString(){
            return Arrays.toString(data);
        }
    }

    /**
     * Your NumArray object will be instantiated and called as such: NumArray obj =
     * new NumArray(nums); obj.update(i,val); int param_2 = obj.sumRange(i,j);
     */
}
