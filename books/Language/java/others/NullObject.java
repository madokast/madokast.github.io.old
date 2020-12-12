package books.language.java.others;

class NullObject {
    public static void main(String[] args) {
        NullObject nullObject = null;
        // false
        System.out.println(nullObject instanceof NullObject);

        String string = null;
        System.out.println(string instanceof String == nullObject instanceof NullObject);

        System.out.println(string==null);

        // boolean b = string==nullObject;

        // System.out.println(b);
    }
}