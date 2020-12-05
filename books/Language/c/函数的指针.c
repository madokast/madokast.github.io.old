#include <stdio.h>

int add(int a, int b) {
    return a + b;
}

int sub(int a, int b) {
    return a - b;
}

int compute(int a, int b, int (*call)(int, int)) {
    return (*call)(a, b);
}

int main() {
    int (*fun)(int, int) = add;

    printf("%d\n", (*fun)(1, 2));

    printf("%d\n",compute(2,3,add));
    printf("%d\n",compute(10,20,sub));

    return 0;
}