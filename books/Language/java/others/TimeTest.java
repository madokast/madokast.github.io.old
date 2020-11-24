package books.language.java.others;

public class TimeTest {
    public static void main(String[] args) {
        for (int i = 0; i < 1000; i++) {
            long s = System.currentTimeMillis();
            test();
            long e = System.currentTimeMillis();
            System.out.println(e-s);
        }
    }

    static void test(){
        int max = 1000*1000;
        int i;
        float temp = 0.0f;
        for(i=0;i<max;i++){
            temp += Math.sin((float)i);
        }
        System.out.println(temp);
    }
}
