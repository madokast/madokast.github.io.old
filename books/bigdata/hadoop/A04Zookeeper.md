# Zookeeper

## 概述

- 定义：分布式协调服务框架，用于解决分布式系统一致性问题

- 本质：分布式文件系统，可以理解为一个数据库

- 应用场景

数据订阅/发布：推拉均可，节点可以被关注，节点改变时，客户可以收到通知。

命名服务

分布式协调通知：心跳检测、工作进度汇报、系统调度。（都是通过节点完成的）

分布式锁：写锁 /exclusive_lock/lock、读锁 /shared_lock/lock

分布式队列：任务按照顺序执行


## Znode：树型层级结构，即是文件夹，也可以存储数据

## 架构：主从架构

- 3 种角色 leader follower observer

<img src="./img/A04zookeeper架构.jpg"></img>

可以看到，从主机可以处理 read，但是写请求都会转发给 leader

<img src="./img/A04zookeeper角色.jpg"></img>

## leader 选举机制

- 时机：服务器启动、leader 挂了

# 安装配置 zookeeper


