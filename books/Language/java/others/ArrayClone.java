package books.language.java.others;

import java.util.Arrays;

public class ArrayClone {
    public static void main(String[] args) {
        int[] arr = {1,2,3};
        System.out.println(Arrays.toString(arr));
        int[] arrCloned = arr.clone();
        arr[1] = 10;
        System.out.println(Arrays.toString(arr));
        System.out.println(Arrays.toString(arrCloned));
    }
    
}
