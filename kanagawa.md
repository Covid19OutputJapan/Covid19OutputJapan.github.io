---
title: "Covid-19 and Output in Japan: Kanagawa"
keywords: sample homepage
tags: [kanagawa]
sidebar: home_sidebar
permalink: kanagawa.html
summary:
---

{% assign fig_loc = "./archives/20210126/Figures/Kanagawa/" %}

## Last update on January 26, 2021

Replications files are available [here](https://github.com/Covid19OutputJapan/Covid19OutputJapan.github.io/tree/main/archives/).

<!--
Link to other Kanagawa pages:
<table>
<tr>
{% assign cnt = 0 %}
{% for page1 in site.pages %}
    {% for tag1 in page1.tags %}
        {% if tag1 == "kanagawa" and page1.name != page.name %}
            <td><a href="{{page1.url | remove: "/" }}">{{page1.permalink}}</a></td>
            {% assign cnt = cnt | plus:1 %}
        {% endif %}
    {% endfor %}
{% endfor %}
</tr>
</table>
-->

### 1. Baseline scenario

{: align="center"}
|![Baseline]({{ fig_loc }}BaselineDecline.png)|

Source: Authors’ calculation.<br>
Note:	See Fujii and Nakata (2021) for a detailed discussion of the scenario.
