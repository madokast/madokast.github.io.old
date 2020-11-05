import sys

print(sys.argv)

import os

print(os.getcwd())

A06 = open(r'books\python\beginner\A06文件.py')

print(A06)

print(A06.read())

A06.close()

print("---------")

with open(r'books\python\beginner\A06文件.py') as steam:
    print(steam.read())

print(__file__)