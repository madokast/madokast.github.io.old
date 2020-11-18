# 进程

- Linux 下 python 使用 os.fork() 创建进程

- win 下使用 multiprocessing 模块创建进程。from multiprocessing import Process # Process(target=task2()).start()

- 进程不共享全局变量

- 进程传消息，自定义进程重写run方法

- 进程池 multiprocessing.pool 

- 进程通讯 queue = multiprocessing.Queue(5) 


# 线程

- threading 模块

- 全局解释器锁——目的：全局更新对象的引用计数