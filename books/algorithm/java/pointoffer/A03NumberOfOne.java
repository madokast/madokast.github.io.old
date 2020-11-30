package books.algorithm.java.pointoffer;

/**
 * 二进制数中 1 的个数 两种解法
 * 第一种 numberOfOne_2 很久以前背过的，错位相加法
 * 
 * 第二种 书上写的，减一&原数
 * 一个数减去 1 后与上原数，则相当于原数最后一个 1 变成 0
 */

public class A03NumberOfOne {
    public static void main(String[] args) {
        test();
    }

    public static void test() {
        A03NumberOfOne a = new A03NumberOfOne();

        System.out.println(a.numberOfOne(9));
        System.out.println(a.numberOfOne(0x7fffffff));
        System.out.println(a.numberOfOne(0xffffffff));
        System.out.println(a.numberOfOne(0x80000000));

        System.out.println("--------------");

        System.out.println(a.numberOfOne_2(9));
        System.out.println(a.numberOfOne_2(0x7fffffff));
        System.out.println(a.numberOfOne_2(0xffffffff));
        System.out.println(a.numberOfOne_2(0x80000000));

        long s;
        int r;

        for (int j = 0; j < 3; j++) {
            s = System.currentTimeMillis();
            r = 0;
            for (int i = 0; i < 0x0fffffff; i++) {
                r += a.numberOfOne(i);
            }
            System.out.println(System.currentTimeMillis() - s);
            System.out.println(r);

            s = System.currentTimeMillis();
            r = 0;
            for (int i = 0; i < 0x0fffffff; i++) {
                r += a.numberOfOne_2(i);
            }
            System.out.println(System.currentTimeMillis() - s);
            System.out.println(r);

        }
    }

    /**
     * 用时 1892 ms
     * 
     * @param a
     * @return
     */
    public int numberOfOne(int a) {
        int c = 0;
        while (a != 0) {
            a = (a - 1) & a;
            c++;
        }
        return c;
    }

    /**
     * 用时 550 s
     * 
     * @param a
     * @return
     */
    public int numberOfOne_2(int a) {
        a = (a & 0x55555555) + ((a & 0xaaaaaaaa) >>> 1);
        a = (a & 0x33333333) + ((a & 0xcccccccc) >>> 2);
        a = (a & 0x0f0f0f0f) + ((a & 0xf0f0f0f0) >>> 4);
        a = (a & 0x00ff00ff) + ((a & 0xff00ff00) >>> 8);
        a = (a & 0x0000ffff) + ((a & 0xffff0000) >>> 16);

        return a;
    }

}
