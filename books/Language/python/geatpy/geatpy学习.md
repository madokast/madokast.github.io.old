# geatpy

## 简介

- Geatpy是一个高性能实用型进化算法工具箱。

- 官网：http://geatpy.com/index.php/home/ 

## 基本结构

- 四大类。Problem问题类、Algorithm算法模板类、Population种群类和PsyPopulation多染色体种群类

### Problem 定义了与问题相关的一些信息

Problem类定义了与问题相关的一些信息，如问题名称name、优化目标的维数M、决策变量的个数Dim、决策变量的范围ranges、决策变量的边界borders等。

maxormins是一个记录着各个目标函数是最小化抑或是最大化的Numpy array行向量。

varTypes是一个记录着决策变量类型的行向量，其中的元素为0表示对应的决策变量是连续型变量；为1表示对应的是离散型变量。

待求解的目标函数定义在aimFunc()的函数中。

**在实际使用时，通过定义一个继承Problem的子类来完成对问题的定义。**

### Population 一个表示种群的类

一个种群包含很多个个体，而每个个体都有一条染色体(若要用多染色体，则使用多个种群、并把每个种群对应个体关联起来即可)。

### PsyPopulation 继承了Population的支持多染色体混合编码的种群类

一个种群包含很多个个体，而每个个体都有多条染色体。

### Algorithm 进化算法的核心类

## 快速入门

见 A02_DTLZ1.py





