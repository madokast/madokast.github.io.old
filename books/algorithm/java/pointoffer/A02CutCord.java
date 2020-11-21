package books.algorithm.java.pointoffer;

/**
 * 绳子长度 len 切成 part 份 问切完后绳子每段长度相乘的最大值 max(Πlen[i])
 */

public class A02CutCord {
    public static void main(String[] args) {
        int ans = solve(8, 3);
        System.out.println(ans);
    }

    private static int solve(int length, int partNum) {
        if (partNum == 1)
            return length;
        else {
            int max = 0;
            for (int i = 1; i < length; i++) {
                max = Math.max(max, solve(length - i, partNum - 1) * i);
            }
            return max;
        }
    }
}
