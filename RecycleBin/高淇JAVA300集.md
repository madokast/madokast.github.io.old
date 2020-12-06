高淇JAVA300集
001~009前言
001如何学习JAVA300集
预科14集。
重代码，重底层，重项目。
遇到难点不纠结，一个月后再看难点。以建立知识体系为目的学习。
002计算机发展史及其未来方向
计算机已成为生活中不可或缺的一部分。
算法是计算机的灵魂。
计算机语言的三代发展。
找工作：不用担心开发人才饱和，因为人越多，软件就越复杂，反而需要更多人才。
未来30年是软件人才的世界。
003常见编程语言介绍
C/C++/JAVA	PHP-webObject-c/Swift	JavaScript/H5Node.JS	PythonC#/.NET	FrotranBasic	…
004JAVA的发展历史和未来
从语言变成了生态体系，因此很难被取代。
1991年，SUN，消费类电子产品的编程（中立，独立于操作系统、CPU）
1998年，JAVA1.2 革命性版本。
JAVA的发展和互联网的发展息息相关。
005JAVA的核心优势和生态体系
核心优势：跨平台。
赶上了互联网的大爆发。计算机界的“英语”
006JAVA三版本
EE包含SE，ME是SE一部分再加自己的东西。
JAVASE（standard edition）定位于个人计算机，发展最差。但是学习应从此开始。
JAVAEE（enterprise edition）定位于企业服务器端。包含完整的SE
JAVAME（micro edition）定位于消费型电子产品，趋于消亡。
注意：JAVAME≠Android
007JAVA的特性和优势
跨平台，可移植。e.g. C/CPP种int可能16位，可能32位，而JAVA都是32位。面向对象。简单性。入门简单。高性能。分布式。多线程。健壮性。（自己崩溃不会危及整个系统）
008JAVA应用程序的运行机制
结合了编译、解释两过程的语言。
源文件(.java)——(JAVA编译器)——字节码文件(.class)送到JRE中，(类装载器)(字节码检验器)(解释器JVM)最后运行于系统平台。
009JDK、JRE、JVM作用和区别
是一种包含关系。
JVM——运行字节码的虚拟机。（光有它不够）
JRE=JVM+库函数+运行Jave应用的必要文件。（有JRE才能运行Java程序）
JDK=JRE+编译器+调试器+其他开发工具
英文：
java virtual machine
java running environment
java development kit
010~022第一個程序
010JDK安装
安装过程略。
关于安装后文件目录：
./bin 二进制文件。exe，dll文件。javac.exe java.exe在这里面
./db 数据文件
./include 头文件
./jre
./lib 库 jar包（字节码打包）src.zip jdk库的源码
011JDK环境变量PATH设置 classpath问题
先新建变量名 JAVA_HOME 地址为JDK的安装文件夹然后设置Path变量 
新增地址 %JAVA_HOME%\bin
**理由：**以后安装的JAVA相关软件会搜索JAVA_HOME
变量JDK1.5以后不同设置classpath，因为JRE会自动搜索当前路径下的class文件、jar文件
012控制台测试JDK安装和配置成功
cmd >>java -version
013写出第一个JAVA程序
014第一个JAVA程序错误总结
遇到错误不要怕，这是你成长的机会。
015第一个JAVA程序详细解释
016常见的DOS命令
cd	dir	cls	上下键	Tab
017常用的Java开发工具
Notpad++	UltraEdit	EditPlus
018eclipse开发环境使用原因
无论是notepad还是IDE，关键是体验带编程的乐趣。
019下载eclipse版本选择和使用
注意和JDK配套，安装32位或64位。
workspace 代码等存放位置
020eclipse下建立JAVA项目_项目的结构
021使用eclipse开发和运行JAVA程序
022开发桌球小游戏
023~045基本數據類型/運算符/類型轉化
023注释
//单行注释
/*
多行注释
/
/*
文档注释。包含注释和JavaDoc标签，生成文档API
*/
024标识符规则
开头：字母、下划线、美元符号$。不能使用数字。
类名：首字母大写
方法名、变量名：首字母小写，第二个单词开始大写
025Java关键字和保留字
026变量的本质
变量代表一个“可操作的存储空间”，位置固定，但是里面的值不确定。
变量需要声明，指定空间大小，读取方式。（强类型语言）
声明与初始化。
027变量的分类
局部分类、成员变量、静态变量。
局部变量：在方法和语句块内定义，生命周期从声明到方法/语句块结束。
成员对象：生命周期伴随对象。会有默认初始值。
静态变量：static修饰符。声明周期从类加载到卸载。（疑问：加载卸载是什么？我可以控制吗？）
028常量constant
固定的值，以及final修饰的变量。
029基本数据类型
相对引用数据类型而言。
三大类，八种。
byte	short	int	long
float	double
char	boolean
030整形变量和整形常量
byte	1字节	-128~127
short	2字节	-2^15~2^15-1
int	4字节	21亿左右
long	8字节	
不同进制。0nn八进制	0xnn十六进制	0bnn二进制
整形常量默认是int类型，nnnnnL->long类型
031浮点型变量/常量
float	4字节	±E38左右	精确度大于7位
double	8字节	±E308左右
浮点常量默认是double，nnF->float类型
不精确，不能用于比较==。精确下使用java.math.bigDecimal
032字符型变量/常量
char	2字节	6w多字符
编码表示’\u0061’==‘a’
转义字符’…’注意
字符串和字符的区别。
033boolean类型/常量
不能用01赋值，要用true/false
034运算符介绍
035算数运算符
二元运算符。±*/%
整数运算，有long时，结果为long，否则一律为int类型
浮点运算，都是float时，结果才是float，否则一律为double
自增自减运算符。注意前后缀问题。
注意类型转化的问题。
036赋值运算符和复制扩展运算符
=	+=	*=	…
037关系运算符
==	!=	>	<	>=	<=
注意==和=区别
注意运用的对象，基本数据类型/引用数据类型
038逻辑运算符
&	|	&&	||	!	^
逻辑与 逻辑或 短路与 短路或
039位运算符
~	&	|	^	<<	>>
位反	位与	位非	位异或	移位
040字符串连接符+
041条件运算符x?y:z
042运算符优先级的问题
043自动类型转化
容量小的可以自动转或为容量大的。有的是有精度损失的。
如long可以自动转化为float
特例：使用常量时，byte b = 12; 只要不超出大小可以接受
044强制类型转化
045溢出问题
046~059語句/方法/带标签的break和continue  
046 Scanner获得键盘输入  
import java.util.Scanner;
Scanner sc = new Scanner(System.in);
String str = sc.nextLine();
047 控制语句介绍  
048 控制语句 if单语句结构  
049 if-else双选择结构  
050 if-else多选择结构  
051 switch多选择结构  
swith(var){case var1:{..break;}....default:...}
052 while循环详解  
053 do-while循环 for循环  
054 嵌套循环  
055 break&continue  
056 带标签的break和continue  
outer: for(int i=0;i<5;i++)
{    for(int j=0;j<5;j++)
	{        if(i==2)            continue outer;            //break outer;//结束外层循环，内层当然也结束了            System.out.println(""+i+j);    }}
057 语句块和方法  
058 方法的重载  
059 递归  
060~067面向對象/垃圾回收
060 面向过程和面向对象  
061 对象是什么
062 对象和类的关系
063 UML入门
064 内存分析详解  
JAVA虚拟机的内存分为三个区域，栈stack，堆heap，方法区method area
栈：①每调用方法就创造一个站帧（存储局部变量、操作数、方法出口）
②每个线程一个栈。线程私有，不可线程间共享。
③先入后出。内存连续。
堆：①存放建好的对象和数组。
②只有一个堆，线程共享。
③不连续。
方法区，又称静态区：①只有一个方法区，线程共享。
②和堆类似，但是只存放类、常量相关的信息，即不变的、唯一的内容。
065 构造方法重载  
066 垃圾回收机制garbage collection GC  
算法：①发现无用的对象。②回收无用对象占领的内存空间。
引用计数法--缺点：循环引用无法识别
引用可达法
067 通用分代垃圾回收详解  
年轻代 年老代 持久代
堆内存分为：Eden Survivor tenured/Old空间
①所有新生成的对象都放在Eden区，这个区目标是快速回收生命周期短的对象，由Minor GC负责。
    操作频繁，效率较高，但是会浪费内存空间。
②Eden区满了，或者N(默认15)次回收后还存在的对象，都放到Old区。
    Old区满了。启动Major GC和Full GC“大扫除”
③持久代对垃圾回收无显著影响。

Minor GC 清理Eden区域，有用对象放到Survivor1或Survivor2区
Major GC 清理Old区。
Full GC 全部清理，对性能有影响。

Survivor1和Survivor2区，存分年轻代。两个区大小相同，每次只有一个在用，另一个为空。

过程：
新对象都在Eden区，满了后，Minor GC清理，还存活的放到Survivor区
Survivor区15次清理还在的对象，放到Old区
Old达到一定比例，执行Major GC；满了后，执行Full GC

JVM调优一般都是针对Full GC，尽量不要调用。

常见内存泄漏：
创建大量无用对象：大量字符串拼接，使用String而非StringBuilder
静态集合类：因为它永远不会被回收，导致集合中的Objec也不会被回收
连接对象未关闭：如IO流 数据库连接对象 网络连接对象
释放对象时没有删除相应的监听器
068~
068 this的本質/對象創建的過程
1.分配空間，初始化
2.顯式初始化
3.構造方法
4.ruturn對象引用

this(...)	調用另一個構造器
this.abc	調用本對象的變量/方法
069 static關鍵字/靜態對象方法
static修飾的方法和對象從屬於類
普通對象和方法從屬於對象
070 靜態初始化塊/繼承樹的向上追溯
static{}
071 Java參數傳值機制
所有Java參數都是傳值（複製一份）
傳引用時，傳地址，雖然是複製的。但是找到的是同一個對象。
072 Java包機制/package
域名倒着写。类名冲突。静态导入。  
073 import詳解
java.util.Data
java.sql.Data
074 繼承/instanceOf
075 方法的重寫override
子類自身的行爲代替父類的行爲
返回值：子類≤父類
	Student Person
訪問權限：子類≥父類
076 Object類/重写toString()
toString()方法
public String toString() {    return getClass().getName() + "@" + Integer.toHexString(hashCode());}
077 equals方法的重写
‘==’->基本类型-值相同，引用类型-地址相同
重写：
1.判断地址this==o
2.判断非空o==null
3.判断类型getClass() !=o.getClass()
4.强制类型转换后比较各属性
078 super父类引用 继承树追溯
super的意义，引用对象中被子类覆盖的属性和方法

构造对象时，总是先调用super()，层层往上一直追溯到Object类
079 封装的使用
private 只有自己类中可以用，子类也不行
default或者不写 同一个包里可以访问（package）