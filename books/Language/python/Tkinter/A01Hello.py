import tkinter
import time
import threading

def set_label(lable,text):
    lable.config(text=text)

def countdown(lable,end_time)->str:
    while True:
        t = time.time()
        set_label(lable,str(end_time-t) + 's')
        time.sleep(0.1)
    


if __name__ == "__main__":
    # 初始化Tk()
    top = tkinter.Tk()
    top.title('hello tkinter')
    l = tkinter.Label(top,text='main')
    lp = l.pack()

    print(l)
    print(lp)

    threading.Thread(target=countdown,args=(l,30+time.time(),)).start()

    # 进入消息循环
    top.mainloop()