routers = {} # 全局变量 注册的路由-方法
def router(path):
    def wrap(method):
        routers[path] = method  # 注册路由
        def core(*args, **kwargs):
            print(f"访问路径{path}，携带参数{args}和{kwargs}")  # 方法增强
            method(*args, **kwargs)  # 核心业务
        return core
    return wrap

@router('/node1')
def n1(name):
    print(f'核心业务，处理节点1，来人是{name}')

@router('/node2')
def n2():
    print(f'核心业务，处理节点2')

# 网关
def gateway(path,*args,**kwargs):
    m =  routers.get(path)
    if m is None:
        print(f'错误，没有可以处理{path}的方法')
    else:
        m(*args,**kwargs)


if __name__ == "__main__":
    # 看看路由注册情况
    print(routers)
    
    # 外部网络访问
    gateway('/node1','zrx')
    gateway('/node1','mdk')
    gateway('/node2')
    gateway('/node3')