package books.algorithm.java.leetcode.L300;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * 301. 删除无效的括号 删除最小数量的无效括号，使得输入的字符串有效，返回所有可能的结果。
 * 
 * 说明: 输入可能包含了除 ( 和 ) 以外的字符。
 * 
 * 示例 1:
 * 
 * 输入: "()())()" 输出: ["()()()", "(())()"] 示例 2:
 * 
 * 输入: "(a)())()" 输出: ["(a)()()", "(a())()"] 示例 3:
 * 
 * 输入: ")(" 输出: [""]
 */

public class L0301RremoveInvalidParentheses {
    public static void main(String[] args) {
        var l = new L0301RremoveInvalidParentheses();
        l.run();
    }

    public void run() {
        System.out.println(removeInvalidParentheses("()())()")); // ["()()()", "(())()"]
        System.out.println(removeInvalidParentheses("(a)())()")); // ["(a)()()", "(a())()"]
        System.out.println(removeInvalidParentheses(")(")); // []
    }

    public List<String> removeInvalidParentheses(String s) {
        int rights = 0;
        int lefts = 0;
        char[] sc = s.toCharArray();
        for (char c : sc) {
            if (c == '(')
                lefts++;
            else if (c == ')') {
                if (lefts == 0)
                    rights++;
                else
                    lefts--;
            }
        }

        System.out.println("需要删除个" + lefts + "左括号，" + rights + "个右括号");

        int len = sc.length;
        boolean[] removed = new boolean[len];
        Set<String> ret = new HashSet<>();

        dps(sc, removed, 0, ret, lefts, rights);

        return new ArrayList<>(ret);
    }

    private void dps(char[] sc, boolean[] removed, int i, Set<String> ret, int leftLeft, int leftRight) {
        if (validate(sc, removed)) {
            ret.add(make(sc, removed));
        } else if ((leftLeft > 0 || leftRight > 0) && i < sc.length) {
            if (sc[i] == '(' && leftLeft > 0) {
                removed[i] = true;
                dps(sc, removed, i + 1, ret, leftLeft - 1, leftRight);
                removed[i] = false;
            } else if (sc[i] == ')' && leftRight > 0) {
                removed[i] = true;
                dps(sc, removed, i + 1, ret, leftLeft, leftRight - 1);
                removed[i] = false;
            }

            dps(sc, removed, i + 1, ret, leftLeft, leftRight);
        }
    }

    private String make(char[] sc, boolean[] removed) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < sc.length; i++) {
            if (!removed[i]) {
                sb.append(sc[i]);
            }
        }

        return sb.toString();
    }

    private boolean validate(char[] sc, boolean[] removed) {
        int lefts = 0;
        for (int i = 0; i < sc.length; i++) {
            if (!removed[i]) {
                if (sc[i] == '(')
                    lefts++;
                else if (sc[i] == ')') {
                    lefts--;
                    if (lefts < 0)
                        return false;
                }
            }
        }

        return lefts == 0;
    }
}
