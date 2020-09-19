## 术语

- **标签 tag**。尖括号+关键字。如&lt;h1&gt;
- **元素**。开始标签+内容+结束标签。如&lt;h1&gt;Java&lt;/h1&gt;
- **空元素**。内容为空的元素。如&lt;br&gt;&lt;/br&gt;
- **自闭合标签**。空元素在开始标签中进行关闭。&lt;br /&gt;
- **虚元素**。不允许假如内容的元素。如&lt;br&gt;
- **属性**。提供了有关 HTML 元素的更多的信息。属性在元素的开始标签中规定，格式为name="value"
- **布尔属性**。使用布尔属性时，只需要写属性名。如&lt;h1 hidden &gt;Java&lt;/h1&gt;
- **自定义属性**。HTML5 规定自定义属性为`data-*`，为了避免和未来新增的属性冲突。

`关键字区分大小写吗？——不区分`

## HTML 文档

- &lt;!DOCTYPE html&gt; 告诉浏览器这是一个 HTML 文档，其中的 html 为布尔属性。
- **元数据**。&lt;head&gt;中的内容。如&lt;title&gt;
- **内容**。&lt;body&gt;中的内容。
- 元素关系。父子元素、后代元素、兄弟元素。
- **HTML 实体**。HTML 中有些字符具有特殊含义（如&lt;&gt;），如果要作为文本输出，应使用实体名称（&amp;lt;&amp;gt;），或编号(&amp;#60;&amp;#62;)。


## 全局属性

- accesskey。定义元素快捷键，用户可以通过Tab+key快捷键快速访问。
- class。
- contenteditable。元素内容是否可编辑。
- contextmenu。元素的上下文菜单。上下文菜单在用户点击元素时显示。只有 Firefox 支持。
- dir。文字方向。
- draggable / dropzone。元素是否可拖动、被拖动数据时是否进行复制、移动或链接。
- hidden。隐藏元素。（**hidden == display:none ≠ visibility:hidden**）
- id。
- lang。元素内容的语言。
- spellcheck。是否对元素进行拼写和语法检查。
- style。行内样式。
- tabindex。Tab次序。
- title。额外信息，鼠标悬浮显示提示。
- translate。无浏览器支持。

# 标签/元素

## 文档和元数据

DOCTYPR、html、head、body、link、meta、script、noscript、style、title

## 文本元素
a、b、br、code、del、em、i、small、sub、sup

## 分组元素
blockquote、div、li、ol、ul、p、pre

## 内容元素（主要HTML5）
h1~h6、article、aside、footer、header、section

## 表格元素
略

## 表单元素
略

## 嵌入元素
audio、canvas、img、video