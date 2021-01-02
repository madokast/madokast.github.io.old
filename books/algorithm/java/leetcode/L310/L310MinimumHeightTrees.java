package books.algorithm.java.leetcode.L310;

import java.util.ArrayList;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * 树是一个无向图，其中任何两个顶点只通过一条路径连接。 换句话说，一个任何没有简单环路的连通图都是一棵树。
 * 
 * 给你一棵包含 n 个节点的数，标记为 0 到 n - 1 。给定数字 n 和一个有 n - 1 条无向边的 edges 列表（每一个边都是一对标签），其中
 * edges[i] = [ai, bi] 表示树中节点 ai 和 bi 之间存在一条无向边。
 * 
 * 可选择树中任何一个节点作为根。当选择节点 x 作为根节点时，设结果树的高度为 h 。在所有可能的树中，具有最小高度的树（即，min(h)）被称为
 * 最小高度树 。
 * 
 * 请你找到所有的 最小高度树 并按 任意顺序 返回它们的根节点标签列表。
 * 
 * 树的 高度 是指根节点和叶子节点之间最长向下路径上边的数量。
 * 
 * 来源：力扣（LeetCode） 链接：https://leetcode-cn.com/problems/minimum-height-trees
 * 著作权归领扣网络所有。商业转载请联系官方授权，非商业转载请注明出处。
 */

public class L310MinimumHeightTrees {
    public static void main(String[] args) {
        L310MinimumHeightTrees l = new L310MinimumHeightTrees();
        System.out.println(l.findMinHeightTrees(4, new int[][]{
            new int[]{1,0},
            new int[]{1,2},
            new int[]{1,3}
        }));
    }

    public List<Integer> findMinHeightTrees(int n, int[][] edges) {
        if(n==1) return List.of(0);
        List<Integer> ret = new ArrayList<>();
        int minHeight = n;
        Set<Integer> nodeSet = new HashSet<>();
        for (int i = 0; i < n; i++) {
            nodeSet.add(i);
        }
        Map<Integer, List<Integer>> paths = new HashMap<>();
        for (int[] path : edges) {
            paths.putIfAbsent(path[0], new ArrayList<>());
            paths.putIfAbsent(path[1], new ArrayList<>());
            paths.get(path[0]).add(path[1]);
            paths.get(path[1]).add(path[0]);
        }

        for (int root = 0; root < n; root++) {
            int maxHeight = bfs(root, paths);
            if (maxHeight < minHeight) {
                ret.clear();
                ret.add(root);
                minHeight = maxHeight;
            } else if (maxHeight == minHeight) {
                ret.add(root);
            }
        }

        return ret;
    }

    private int bfs(int root, Map<Integer, List<Integer>> paths) {
        int height = 0;
        Set<Integer> visited = new HashSet<>();
        visited.add(root);

        Deque<Integer> queue = new LinkedList<>();
        for (Integer c : paths.get(root)) {
            queue.push(c);
        }

        while (!queue.isEmpty()) {
            Set<Integer> nextQueue = new HashSet<>();
            while (!queue.isEmpty()) {
                Integer c = queue.pop();
                visited.add(c);
                for (Integer cc : paths.get(c)) {
                    if (!visited.contains(cc))
                        nextQueue.add(cc);
                }
            }

            queue = new LinkedList<>(nextQueue);
            height++;
        }

        return height;
    }
}
