var MarkdownRender = {
    CONSTENT: {
        HTML_H1: ['<h1>', '</h1>'],
        HTML_H2: ['<h2>', '</h2>'],
        HTML_H3: ['<h3>', '</h3>'],
        HTML_P: ['<p>', '</p>'],
        HTML_PRE: ['<pre>', '</pre>'],
        HTML_UL: ['<ul>', '</ul>'],
        HTML_OL: ['<ol>', '</ol>'],
        HTML_LI: ['<li>', '</li>'],
        HTML_STRONG: ['<strong>', '</strong>'],
        HTML_EM: ['<em>', '</em>'],
        HTML_INLINE_CODE: ['<span class="markdown-inline-code" >', '</span>'],
        MD_STRONG: '**',
        MD_EM: '*',
        MD_INLINE_CODE: '`',
        MARKDOWM_TITLE: 1,
        MARKDOWM_NORMAL: 2,
        MARKDOWM_EMPTY_LINE: 3,
        MARKDOWM_CODE: 4,
        MARKDOWM_ORDER_LIST: 5,
        MARKDOWM_UNORDER_LIST: 6,
        MARKDOWM_HTML: 7,
        HTMLS: undefined // json 内无法引用内部变量 https://www.imooc.com/wenda/detail/571609 
    },
    markDownUrl: '',
    baseUrl: '',
    // 发送请求，渲染 markdown
    // url markdown文档地址
    // container 渲染目的地dom
    render: function (url, container) {
        this.markDownUrl = url; // https://madokast.github.io/draft/markdown-learn.md
        this.baseUrl = url.substring(0, url.lastIndexOf('/'));
        var that = this;

        //步骤一:创建异步对象
        var ajax = new XMLHttpRequest();
        //步骤二:设置请求的url参数,参数一是请求的类型,参数二是请求的url,可以带参数,动态的传递参数starName到服务端
        ajax.open('get', url);
        //步骤三:发送请求
        ajax.send();
        //步骤四:注册事件 onreadystatechange 状态改变就会调用
        ajax.onreadystatechange = function () {
            if (ajax.readyState == 4 && ajax.status == 200) {
                //步骤五 如果能够进到这个判断 说明 数据 完美的回来了,并且请求的页面是存在的
                that.doRender(ajax.responseText, container);
            }
        }
    },
    doRender: function (content, container) {
        var lines = content.split('\n');
        var content = [];
        var levelStack = [];
        for (var i = 0; i < lines.length; i++) {
            var line = lines[i].trim();
            levelStack = this.renderLine(line, content, levelStack);
        }
        container.innerHTML = content.join('\n');
    },
    // 渲染一行
    renderLine: function (line, content, levelStack) {
        var stateArr = this.resolve(line);
        var curMarkDown = stateArr[0];

        // 多行结构体的前后标记
        if (curMarkDown == this.CONSTENT.MARKDOWM_CODE || curMarkDown == this.CONSTENT.MARKDOWM_ORDER_LIST
            || curMarkDown == this.CONSTENT.MARKDOWM_UNORDER_LIST) {
            if (levelStack.length == 0) {
                levelStack[levelStack.length] = stateArr; // push
                content[content.length] = this.CONSTENT.HTMLS[curMarkDown][0];
            } else {
                var peek = levelStack[levelStack.length - 1][0]; // peek
                // MARKDOWM_CODE 需要单独考虑，因为在块结果中，MARKDOWM_ORDER_LIST 和 MARKDOWM_UNORDER_LIST 每行具有标记
                // 但是 MARKDOWM_CODE 只有开始、结束标记
                if (peek == this.CONSTENT.MARKDOWM_CODE && curMarkDown == this.CONSTENT.MARKDOWM_CODE) {
                    levelStack.length = levelStack.length - 1; // pop
                    content[content.length] = this.CONSTENT.HTMLS[curMarkDown][1];
                } else if (peek != this.CONSTENT.MARKDOWM_CODE && curMarkDown != this.CONSTENT.MARKDOWM_CODE && peek != curMarkDown) {
                    levelStack[levelStack.length] = stateArr; // push
                    content[content.length] = this.CONSTENT.HTMLS[curMarkDown][0];
                }
            }
        } else {
            if (levelStack.length > 0) {
                var peek = levelStack[levelStack.length - 1][0]; // peek
                if (peek == this.CONSTENT.MARKDOWM_CODE) {
                    // do nothing
                } else {
                    levelStack.length = levelStack.length - 1; // pop
                    content[content.length] = this.CONSTENT.HTMLS[peek][1];
                }


            }
        }

        content[content.length] = this.solve(line, stateArr, levelStack);
        //console.log(line + " --> " + content[content.length - 1]);

        return levelStack;
    },
    // 解析一行，识别是 MARKDOWM_TITLE 还是 MARKDOWM_NORMAL 等等
    // 输出数组，[0] 表示改行识别码，如 MARKDOWM_TITLE，[1]根据识别码不同含义不同
    resolve: function (line) {
        if (line.length == 0) {
            // 空行
            return [this.CONSTENT.MARKDOWM_EMPTY_LINE];
        } else {
            var firstCh = line[0];
            if (firstCh == '#') {
                // 可能是标题
                var i = 1;
                while (line[i] == '#') {
                    i++;
                }
                if (line[i] == ' ') {
                    return [this.CONSTENT.MARKDOWM_TITLE, i];
                } else {
                    return [this.CONSTENT.MARKDOWM_NORMAL];
                }
            } else if (firstCh == '`') {
                // 可能是代码块
                if (line.length >= 3 && line[1] == '`' && line[2] == '`') {
                    return [this.CONSTENT.MARKDOWM_CODE, line.substring(3)];
                } else {
                    return [this.CONSTENT.MARKDOWM_NORMAL];
                }
            } else if (firstCh == '-') {
                if (line.length >= 3 && line[1] == ' ') {
                    return [this.CONSTENT.MARKDOWM_UNORDER_LIST];
                } else {
                    return [this.CONSTENT.MARKDOWM_NORMAL];
                }
            } else if (!isNaN(firstCh)) {
                // 可能是有序列表
                var index = line.indexOf('.');
                if (index == -1) {
                    return [this.CONSTENT.MARKDOWM_NORMAL];
                } else {
                    var i = 1;
                    while (i < index) {
                        if (!isNaN(line[i])) {
                            i++;
                        } else {
                            break;
                        }
                    }
                    if (i == index) {
                        return [this.CONSTENT.MARKDOWM_ORDER_LIST, parseInt(line.substring(0, index))];
                    } else {
                        return [this.CONSTENT.MARKDOWM_NORMAL];
                    }
                }
            } else if (firstCh == '<') {
                var fisrtGt = line.indexOf('>');
                if (fisrtGt == -1) return [this.CONSTENT.MARKDOWM_NORMAL];
                var secondGt = line.indexOf('>', fisrtGt);
                if (secondGt == -1) return [this.CONSTENT.MARKDOWM_NORMAL];
                else return [this.CONSTENT.MARKDOWM_HTML];
            } else {
                return [this.CONSTENT.MARKDOWM_NORMAL];
            }
        }
    },
    // 渲染一行
    solve: function (line, stateArr, levelStack) {
        var curState = stateArr[0];
        if (curState == this.CONSTENT.MARKDOWM_TITLE) {
            var level = stateArr[1];
            if (level == 1) {
                return this.CONSTENT.HTML_H1[0] + line.substring(level + 1) + this.CONSTENT.HTML_H1[1];
            } else if (level == 2) {
                return this.CONSTENT.HTML_H2[0] + line.substring(level + 1) + this.CONSTENT.HTML_H2[1];
            } else if (level > 2) {
                return this.CONSTENT.HTML_H3[0] + line.substring(level + 1) + this.CONSTENT.HTML_H3[1];
            }
        } else if (curState == this.CONSTENT.MARKDOWM_NORMAL) {
            if (levelStack.length > 0) {
                return line;
            } else {
                return this.CONSTENT.HTML_P[0] + this.solveInLine(line) + this.CONSTENT.HTML_P[1];
            }
        } else if (curState == this.CONSTENT.MARKDOWM_EMPTY_LINE) {
            return '';
        } else if (curState == this.CONSTENT.MARKDOWM_CODE) {
            // do nothing
        } else if (curState == this.CONSTENT.MARKDOWM_ORDER_LIST) {
            var number = stateArr[1];
            return this.CONSTENT.HTML_LI[0] + this.solveInLine(line.substring(("" + number).length + 2)) + this.CONSTENT.HTML_LI[1];
        } else if (curState == this.CONSTENT.MARKDOWM_UNORDER_LIST) {
            return this.CONSTENT.HTML_LI[0] + this.solveInLine(line.substring(2)) + this.CONSTENT.HTML_LI[1];
        } else if (curState == this.CONSTENT.MARKDOWM_HTML) {
            return this.CONSTENT.HTML_P[0] + this.solveHTML(line) + this.CONSTENT.HTML_P[1];
        }
    },
    // 行内渲染
    solveInLine: function (line) {
        var ret = '';
        var boldStart = line.indexOf(this.CONSTENT.MD_STRONG);
        var start = 0;
        while (boldStart != -1) {
            var boldEnd = line.indexOf(this.CONSTENT.MD_STRONG, boldStart + this.CONSTENT.MD_STRONG.length);
            if (boldEnd != -1) {
                ret = ret + line.substring(start, boldStart) + this.CONSTENT.HTML_STRONG[0]
                    + line.substring(boldStart + this.CONSTENT.MD_STRONG.length, boldEnd) + this.CONSTENT.HTML_STRONG[1];
                start = boldEnd + 2;
                boldStart = line.indexOf(this.CONSTENT.MD_STRONG, start);
            } else {
                break;
            }
        }

        line = ret + line.substring(start);
        ret = '';
        start = 0;

        var emStart = line.indexOf(this.CONSTENT.MD_EM);
        while (emStart != -1) {
            var emEnd = line.indexOf(this.CONSTENT.MD_EM, emStart + this.CONSTENT.MD_EM.length);
            if (emEnd != -1) {
                ret = ret + line.substring(start, emStart) + this.CONSTENT.HTML_EM[0]
                    + line.substring(emStart + this.CONSTENT.MD_EM.length, emEnd) + this.CONSTENT.HTML_EM[1];
                start = emEnd + this.CONSTENT.MD_EM.length;
                emStart = line.indexOf(this.CONSTENT.MD_EM, start);
            } else {
                break;
            }
        }

        line = ret + line.substring(start);
        ret = '';
        start = 0;

        var inLineCodeStart = line.indexOf(this.CONSTENT.MD_INLINE_CODE);
        while (inLineCodeStart != -1) {
            var inLineCodeEnd = line.indexOf(this.CONSTENT.MD_INLINE_CODE, inLineCodeStart + this.CONSTENT.MD_INLINE_CODE.length);
            if (inLineCodeEnd != -1) {
                ret = ret + line.substring(start, inLineCodeStart) + this.CONSTENT.HTML_INLINE_CODE[0]
                    + line.substring(inLineCodeStart + this.CONSTENT.MD_INLINE_CODE.length, inLineCodeEnd) + this.CONSTENT.HTML_INLINE_CODE[1];
                start = inLineCodeEnd + this.CONSTENT.MD_INLINE_CODE.length;
                inLineCodeStart = line.indexOf(this.CONSTENT.MD_INLINE_CODE, start);
            } else {
                break;
            }
        }

        line = ret + line.substring(start);

        return line;

    },
    // 渲染 html
    // 主要是渲染 src
    solveHTML: function (line) {
        // <img src="./small.png" alt="图片名称"/>
        var srcIndex = line.indexOf('src');
        if (srcIndex != -1) {
            var srcRef1 = line.indexOf('"', srcIndex);
            var srcRef2 = srcRef1 == -1 ? -1 : line.indexOf('"', srcRef1 + 1);
            if (srcRef2 != -1) {
                var src = line.substring(srcRef1 + 1, srcRef2).trim();
                return line.substring(0, srcRef1 + 1) + this.baseUrl + '/' + src + line.substring(srcRef2);
            }
        }
    }
};

MarkdownRender.CONSTENT.HTMLS = [0, 1, 2, 3, MarkdownRender.CONSTENT.HTML_PRE, MarkdownRender.CONSTENT.HTML_OL, MarkdownRender.CONSTENT.HTML_UL, 7];