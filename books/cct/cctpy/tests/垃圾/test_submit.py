print(__name__)

try:
    from books.cct.cctpy.cctpy import *
except ModuleNotFoundError:
    pass

import sys
sys.path.append(r'C:\Users\madoka_9900\Documents\github\madokast.github.io\books\cct\cctpy')

from cctpy import *

import multiprocessing

def add(x,y):
    print(__name__)
    return x+y


# if __name__ == "__main__":

r = BaseUtils.submit_process_task(
    task=add,
    param_list=[
        [1,2],[3,4]
    ],
    concurrency_level=1
)

print(r)