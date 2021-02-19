打开一个完结的番，犬夜叉 `https://www.bilibili.com/bangumi/play/ss28324/`

点击追番，依次发送两个请求

```
curl 'https://data.bilibili.com/log/web?0000171612493840526https%3A%2F%2Fwww.bilibili.com%2Fbangumi%2Fplay%2Fss28324%2F|666.25.selfDef.click_follow||1612493840000|0|0|498x645|1|\{%22event%22:%22click_follow%22,%22value%22:\{%22season_id%22:28324,%22ep_id%22:289986,%22businesstype%22:1,%22title%22:%22%E7%8A%AC%E5%A4%9C%E5%8F%89%22\},%22bsource%22:%22search_baidu%22\}|%22\{\%22hitGroup\%22:\%22658\%22\}%22|https%3A%2F%2Fspace.bilibili.com%2F|319F1BC7-6501-1BE8-F5BD-04F60E0719A691223infoc|en-US|null' \
  -X 'POST' \
  -H 'authority: data.bilibili.com' \
  -H 'content-length: 0' \
  -H 'pragma: no-cache' \
  -H 'cache-control: no-cache' \
  -H 'sec-ch-ua: "Google Chrome";v="87", " Not;A Brand";v="99", "Chromium";v="87"' \
  -H 'sec-ch-ua-mobile: ?0' \
  -H 'user-agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.141 Safari/537.36' \
  -H 'content-type: text/plain;charset=UTF-8' \
  -H 'accept: */*' \
  -H 'origin: https://www.bilibili.com' \
  -H 'sec-fetch-site: same-site' \
  -H 'sec-fetch-mode: no-cors' \
  -H 'sec-fetch-dest: empty' \
  -H 'referer: https://www.bilibili.com/bangumi/play/ss28324/' \
  -H 'accept-language: en-US,en;q=0.9,zh-CN;q=0.8,zh;q=0.7,ja;q=0.6,zh-TW;q=0.5' \
  -H $'cookie: ■■■ \
  --compressed
```

```
curl 'https://api.bilibili.com/pgc/web/follow/add' \
  -H 'authority: api.bilibili.com' \
  -H 'pragma: no-cache' \
  -H 'cache-control: no-cache' \
  -H 'sec-ch-ua: "Google Chrome";v="87", " Not;A Brand";v="99", "Chromium";v="87"' \
  -H 'accept: application/json, text/plain, */*' \
  -H 'sec-ch-ua-mobile: ?0' \
  -H 'user-agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.141 Safari/537.36' \
  -H 'content-type: application/x-www-form-urlencoded' \
  -H 'origin: https://www.bilibili.com' \
  -H 'sec-fetch-site: same-site' \
  -H 'sec-fetch-mode: cors' \
  -H 'sec-fetch-dest: empty' \
  -H 'referer: https://www.bilibili.com/bangumi/play/ss28324/' \
  -H 'accept-language: en-US,en;q=0.9,zh-CN;q=0.8,zh;q=0.7,ja;q=0.6,zh-TW;q=0.5' \
  -H $'cookie: ■■■' \
  --data-raw 'season_id=28324&csrf=e464a75ec86c0ae02bb9edc128907936' \
  --compressed
```

现在已经知道，同一个番，追番后推掉，再追，发送的请求是一样的，只有第一个请求，url中两处时间戳不同

0000171612493840526

这个数字17后面的是时间戳

1612493840000

这个也是时间戳。


再打开另一个番，超炮，`https://www.bilibili.com/bangumi/play/ss425/`

追番时发送的两个请求是

```
curl 'https://data.bilibili.com/log/web?0000171612494029710https%3A%2F%2Fwww.bilibili.com%2Fbangumi%2Fplay%2Fss425%2F|666.25.selfDef.click_follow||1612494029000|0|0|498x645|1|\{%22event%22:%22click_follow%22,%22value%22:\{%22season_id%22:425,%22ep_id%22:84363,%22businesstype%22:1,%22title%22:%22%E6%9F%90%E7%A7%91%E5%AD%A6%E7%9A%84%E8%B6%85%E7%94%B5%E7%A3%81%E7%82%AE%22\},%22bsource%22:%22search_baidu%22\}|%22\{\%22hitGroup\%22:\%22658\%22\}%22|https%3A%2F%2Fspace.bilibili.com%2F|319F1BC7-6501-1BE8-F5BD-04F60E0719A691223infoc|en-US|null' \
  -X 'POST' \
  -H 'authority: data.bilibili.com' \
  -H 'content-length: 0' \
  -H 'pragma: no-cache' \
  -H 'cache-control: no-cache' \
  -H 'sec-ch-ua: "Google Chrome";v="87", " Not;A Brand";v="99", "Chromium";v="87"' \
  -H 'sec-ch-ua-mobile: ?0' \
  -H 'user-agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.141 Safari/537.36' \
  -H 'content-type: text/plain;charset=UTF-8' \
  -H 'accept: */*' \
  -H 'origin: https://www.bilibili.com' \
  -H 'sec-fetch-site: same-site' \
  -H 'sec-fetch-mode: no-cors' \
  -H 'sec-fetch-dest: empty' \
  -H 'referer: https://www.bilibili.com/bangumi/play/ss425/' \
  -H 'accept-language: en-US,en;q=0.9,zh-CN;q=0.8,zh;q=0.7,ja;q=0.6,zh-TW;q=0.5' \
  -H $'cookie: ■■■' \
  --compressed
```

```
curl 'https://api.bilibili.com/pgc/web/follow/add' \
  -H 'authority: api.bilibili.com' \
  -H 'pragma: no-cache' \
  -H 'cache-control: no-cache' \
  -H 'sec-ch-ua: "Google Chrome";v="87", " Not;A Brand";v="99", "Chromium";v="87"' \
  -H 'accept: application/json, text/plain, */*' \
  -H 'sec-ch-ua-mobile: ?0' \
  -H 'user-agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.141 Safari/537.36' \
  -H 'content-type: application/x-www-form-urlencoded' \
  -H 'origin: https://www.bilibili.com' \
  -H 'sec-fetch-site: same-site' \
  -H 'sec-fetch-mode: cors' \
  -H 'sec-fetch-dest: empty' \
  -H 'referer: https://www.bilibili.com/bangumi/play/ss425/' \
  -H 'accept-language: en-US,en;q=0.9,zh-CN;q=0.8,zh;q=0.7,ja;q=0.6,zh-TW;q=0.5' \
  -H $'cookie: ■■■' \
  --data-raw 'season_id=425&csrf=e464a75ec86c0ae02bb9edc128907936' \
  --compressed
```

##  对比第一个 POST

1. refer 头不同

犬夜叉：`https://www.bilibili.com/bangumi/play/ss28324/`

超炮：`https://www.bilibili.com/bangumi/play/ss425/`

不同之处就是当前网页 url。`https://www.bilibili.com/bangumi/play/ss■■■/`

2. 请求 url 不同，我们详细对比一下

犬夜叉：

```
https://data.bilibili.com/log/web?0000171612493840526https%3A%2F%2Fwww.bilibili.com%2Fbangumi%2Fplay%2Fss28324%2F|666.25.selfDef.click_follow||1612493840000|0|0|498x645|1|\{%22event%22:%22click_follow%22,%22value%22:\{%22season_id%22:28324,%22ep_id%22:289986,%22businesstype%22:1,%22title%22:%22%E7%8A%AC%E5%A4%9C%E5%8F%89%22\},%22bsource%22:%22search_baidu%22\}|%22\{\%22hitGroup\%22:\%22658\%22\}%22|https%3A%2F%2Fspace.bilibili.com%2F|319F1BC7-6501-1BE8-F5BD-04F60E0719A691223infoc|en-US|null
```

超炮：

```
https://data.bilibili.com/log/web?0000171612494029710https%3A%2F%2Fwww.bilibili.com%2Fbangumi%2Fplay%2Fss425%2F|666.25.selfDef.click_follow||1612494029000|0|0|498x645|1|\{%22event%22:%22click_follow%22,%22value%22:\{%22season_id%22:425,%22ep_id%22:84363,%22businesstype%22:1,%22title%22:%22%E6%9F%90%E7%A7%91%E5%AD%A6%E7%9A%84%E8%B6%85%E7%94%B5%E7%A3%81%E7%82%AE%22\},%22bsource%22:%22search_baidu%22\}|%22\{\%22hitGroup\%22:\%22658\%22\}%22|https%3A%2F%2Fspace.bilibili.com%2F|319F1BC7-6501-1BE8-F5BD-04F60E0719A691223infoc|en-US|null
```

关注时间戳和还有

season_id 28324 / 425

ep_id 289986 / 84363

## 对比第二个 POST

1. refer 头不同，情况和第一个 POST 一致

2. 发送的 data 不同

犬夜叉：season_id=28324

超炮：season_id=425

## 如果我们只发送第二个请求

答：可以追番成功

## 如果写一个错误的ep_id，或者不发 ep_id 会怎么样？

答：不需要发送第一个请求。所以没有尝试。

## 现在我需要知道一个番的 season_id 和 ep_id

缘之空：`https://www.bilibili.com/bangumi/play/ss36377/`

日在校园：`https://www.bilibili.com/bangumi/play/ss35583/`

