package books.algorithm.java.pointoffer;

/**
 * 类似问题，青蛙跳台阶
 */

public class A01Fibonacci {
    public static void main(String[] args) {
        int r = recur(10);
        int l = loop(10);
        System.out.println(r);
        System.out.println(l);
    }

    private static int recur(int n){
        if(n<=2) return 1;
        else return recur(n-1)+recur(n-2);
    }

    private static int loop(int n){
        if(n<=2) return 1;
        else{
            int pre = 1, cur = 1, nex;
            for (int i = 2; i < n; i++) {
                nex = pre + cur;
                pre = cur;
                cur = nex;
            }

            return cur;
        }
    }
}