package books.algorithm.java.pointoffer;

public class Hello {
    public static void main(String[] args) {
        var l = System.currentTimeMillis();
        var a = f(40);
        var l2 = System.currentTimeMillis();
        System.out.println((l2-l) + "  " + a);
    }

    static int f(int i) {
        if (i < 2)
            return 1;
        else
            return f(i - 1) + f(i - 2);
    }
}
