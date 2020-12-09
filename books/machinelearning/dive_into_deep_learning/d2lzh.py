import random
from mxnet import autograd, nd
from IPython import display
from matplotlib import pyplot as plt

def hello():
    print("hello")

def use_svg_display():
    # 用矢量图显示
    display.set_matplotlib_formats('svg')

def set_figsize(figsize=(3.5, 2.5)):
    use_svg_display()
    # 设置图的尺寸
    plt.rcParams['figure.figsize'] = figsize

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

# Fashion-MNIST 数据标签含义对照
def get_fashion_mnist_labels(labels):
    text_labels = ['t-shirt', 'trouser', 'pullover', 'dress', 'coat',
                   'sandal', 'shirt', 'sneaker', 'bag', 'ankle boot']
    return [text_labels[int(i)] for i in labels]

# Fashion-MNIST 绘制图象和标签
def show_fashion_mnist(images, labels):
    use_svg_display()
    # 这里的_表示我们忽略（不使用）的变量
    _, figs = plt.subplots(1, len(images), figsize=(12, 12))
    for f, img, lbl in zip(figs, images, labels):
        f.imshow(img.reshape((28, 28)).asnumpy())
        f.set_title(lbl)
        f.axes.get_xaxis().set_visible(False)
        f.axes.get_yaxis().set_visible(False)