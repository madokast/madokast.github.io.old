import random
from mxnet import autograd, nd

def hello():
    print("hello")


# 遍历数据集并不断读取小批量数据样本
def data_iter(batch_size, features, labels):
    """
    典型使用方法
    batch_size = 10
    for X, y in d2lzh.data_iter(batch_size, features, labels):
        print(X, y)
        break
    """
    num_examples = len(features)
    indices = list(range(num_examples))
    random.shuffle(indices)  # 样本的读取顺序是随机的
    for i in range(0, num_examples, batch_size):
        j = nd.array(indices[i: min(i + batch_size, num_examples)])
        yield features.take(j), labels.take(j)  # take函数根据索引返回对应元素

# 线性回归矢量计算
def linreg(X, w, b):
    return nd.dot(X, w) + b

# 损失函数
def squared_loss(y_hat, y):
    return (y_hat - y.reshape(y_hat.shape)) ** 2 / 2

# 小批量随机批量下降
def sgd(params, lr, batch_size):
    """
    params 要优化的参数 如 w b
    lr 学习速率
    batch_size 一批数据大小
    """
    for param in params:
        param[:] = param - lr * param.grad / batch_size