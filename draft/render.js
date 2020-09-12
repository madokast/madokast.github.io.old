window.onload = function () {
    var container = document.getElementById('markdown-container');
    var url = container.getAttribute('resource');
    render(url, container);
}

// 发送请求，渲染 markdown
// url markdown文档地址
// container 渲染目的地
function render(url, container) {
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
            doRender(ajax.responseText, container);
        }
    }
};

var MARKDOWM_TITLE = '#'; // 标题
var HTML_H1 = ['<h1>', '</h1>']; var HTML_H2 = ['<h2>', '</h2>']; var HTML_H3 = ['<h3>', '</h3>'];
var HTML_P = ['<p>', '</p>'];
var HTML_PRE = ['<pre>', '</pre>']; var HTML_UL = ['<ul>', '</ul>']; var HTML_OL = ['<ol>', '</ol>']; var HTML_LI = ['<li>', '</li>'];
var MARKDOWM_TITLE = 1; var MARKDOWM_NORMAL = 2; var MARKDOWM_EMPTY_LINE = 3; var MARKDOWM_CODE = 4;
var MARKDOWM_ORDER_LIST = 5; var MARKDOWM_UNORDER_LIST = 6; var MARKDOWM_HTML = 7;
var HTMLS = [0, 1, 2, 3, HTML_PRE, HTML_OL, HTML_UL, 7];

function doRender(content, container) {
    var lines = content.split('\n');
    var content = [];
    var levelStack = [];
    for (var i = 0; i < lines.length; i++) {
        var line = lines[i].trim();
        levelStack = renderLine(line, content, levelStack);
    }
    container.innerHTML = content.join('\n');
};

// 渲染一行
function renderLine(line, content, levelStack) {
    var stateArr = resolve(line);
    var curMarkDown = stateArr[0];

    // 多行结构体的前后标记
    // TODO 
    if (curMarkDown == MARKDOWM_CODE || curMarkDown == MARKDOWM_ORDER_LIST || curMarkDown == MARKDOWM_UNORDER_LIST) {
        if (levelStack.length == 0) {
            levelStack[levelStack.length] = stateArr; // push
            content[content.length] = HTMLS[curMarkDown][0];
        } else {
            var peek = levelStack[levelStack.length - 1][0]; // peek

            // MARKDOWM_CODE 需要单独考虑，因为在块结果中，MARKDOWM_ORDER_LIST 和 MARKDOWM_UNORDER_LIST 每行具有标记
            // 但是 MARKDOWM_CODE 只有开始、结束标记
            if (peek == MARKDOWM_CODE && curMarkDown == MARKDOWM_CODE) {
                levelStack.length = levelStack.length - 1; // pop
                content[content.length] = HTMLS[curMarkDown][1];
            } else if (peek != MARKDOWM_CODE && curMarkDown != MARKDOWM_CODE) {
                levelStack[levelStack.length] = stateArr; // push
                content[content.length] = HTMLS[curMarkDown][0];
            }
        }
    } else {
        if (levelStack.length > 0) {
            var peek = levelStack[levelStack.length - 1][0]; // peek
            if (peek == MARKDOWM_CODE) {
                // do nothing
            } else {
                levelStack.length = levelStack.length - 1; // pop
                content[content.length] = HTMLS[peek][1];
            }


        }
    }

    console.log(line + "--" + levelStack);

    content[content.length] = solve(line, stateArr, levelStack);

    return levelStack;
};

// 解析一行，识别是 MARKDOWM_TITLE 还是 MARKDOWM_NORMAL 等等
// 输出数组，[0] 表示改行识别码，如 MARKDOWM_TITLE，[1]根据识别码不同含义不同
function resolve(line) {
    if (line.length == 0) {
        // 空行
        return [MARKDOWM_EMPTY_LINE];
    } else {
        var firstCh = line[0];
        if (firstCh == '#') {
            // 可能是标题
            var i = 1;
            while (line[i] == '#') {
                i++;
            }
            if (line[i] == ' ') {
                return [MARKDOWM_TITLE, i];
            } else {
                return [MARKDOWM_NORMAL];
            }
        } else if (firstCh == '`') {
            // 可能是代码块
            if (line.length >= 3 && line[1] == '`' && line[2] == '`') {
                return [MARKDOWM_CODE, line.substring(3)];
            } else {
                return [MARKDOWM_NORMAL];
            }
        } else if (firstCh == '-') {
            if (line.length >= 3 && line[1] == ' ') {
                return [MARKDOWM_UNORDER_LIST];
            } else {
                return [MARKDOWM_NORMAL];
            }
        } else if (!isNaN(firstCh)) {
            // 可能是有序列表
            var index = line.indexOf('.');
            if (index == -1) {
                return [MARKDOWM_NORMAL];
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
                    return [MARKDOWM_ORDER_LIST, parseInt(line.substring(0, index))];
                } else {
                    return [MARKDOWM_NORMAL];
                }
            }
        } else if (firstCh == '<') {
            var fisrtGt = line.indexOf('>');
            if (fisrtGt == -1) return [MARKDOWM_NORMAL];
            var secondGt = line.indexOf('>', fisrtGt);
            if (secondGt == -1) return [MARKDOWM_NORMAL];
            else return [MARKDOWM_HTML];
        } else {
            return [MARKDOWM_NORMAL];
        }
    }
};

// 渲染一行，
function solve(line, stateArr, levelStack) {
    var curState = stateArr[0];
    if (curState == MARKDOWM_TITLE) {
        var level = stateArr[1];
        if (level == 1) {
            return HTML_H1[0] + line.substring(level + 1) + HTML_H1[1];
        } else if (level == 2) {
            return HTML_H2[0] + line.substring(level + 1) + HTML_H2[1];
        } else if (level > 2) {
            return HTML_H3[0] + line.substring(level + 1) + HTML_H3[1];
        }
    } else if (curState == MARKDOWM_NORMAL) {
        if (levelStack.length > 0) {
            return line;
        } else {
            return HTML_P[0] + line + HTML_P[1];
        }
    } else if (curState == MARKDOWM_EMPTY_LINE) {
        return '';
    } else if (curState == MARKDOWM_CODE) {
        // do nothing
    } else if (curState == MARKDOWM_ORDER_LIST) {
        var number = stateArr[1];
        return HTML_LI[0] + line.substring(("" + number).length + 2) + HTML_LI[1];
    } else if (curState == MARKDOWM_UNORDER_LIST) {
        return HTML_LI[0] + line.substring(2) + HTML_LI[1];
    } else if (curState == MARKDOWM_HTML) {
        return HTML_P[0] + line + HTML_P[1];
    }
}
