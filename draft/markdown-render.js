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

var img_left = '<img';

function doRender(content, container) {
    var lines = content.split('\n');
    var content = [];
    var start = 0;
    for (var i = 0; i < lines.length; i++) {
        var line = lines[i].trim();
        if (line.startsWith(img_left)) {
            content[content.length] = markdown.toHTML(lines.slice(start, i).join('\n'));
            content[content.length] = line;
            start = i + 1;
        }
    }
    content[content.length] = markdown.toHTML(lines.slice(start, i).join('\n'));
    container.innerHTML = content.join('\n');
};


window.onload = function () {
    var container = document.getElementById('container');
    var url = container.getAttribute('url');
    render(url, container);
}
