---
title: "Covid-19 and Output in Japan: Tokyo"
keywords: sample homepage
tags: [tokyo]
sidebar: home_sidebar
permalink: tokyo.html
summary:
---

{% assign fig_loc = "./archives/20210126/Figures/Tokyo/" %}

## Last update on January 26, 2021

Replications files are available [here](https://github.com/Covid19OutputJapan/Covid19OutputJapan.github.io/tree/main/archives/).

Link to other Tokyo pages:
<table>
<tr>
{% assign cnt = 0 %}
{% for page1 in site.pages %}
    {% for tag1 in page1.tags %}
        {% if tag1 == "tokyo" and page1.name != page.name %}
            <td><a href="{{page1.url | remove: "/" }}">{{page1.permalink}}</a></td>
            {% assign cnt = cnt | plus:1 %}
        {% endif %}
<!--
        {% if cnt == 1 %}
            <td>here</td>
            {% assign cnt = 0 %}
        {% endif %}
-->
    {% endfor %}
{% endfor %}
</tr>
</table>

### 1. Baseline scenario

{: align="center"}
|![Baseline]({{ fig_loc }}BaselineDecline.png)|

Source: Authors’ calculation.<br>
Note:	See Fujii and Nakata (2021) for a detailed discussion of the scenario.

<!--
### 2. Rapid-decline scenario

{: align="center"}
|![Rapid]({{ fig_loc }}RapidDecline.png)|

Source: Authors’ calculation.<br>
Note:	See Fujii and Nakata (2021) for a detailed discussion of the scenario.

### 3. Gradual-decline scenario

{: align="center"}
|![Gradual]({{ fig_loc }}GradualDecline.png)|

Source: Authors’ calculation.<br>
Note:	See Fujii and Nakata (2021) for a detailed discussion of the scenario.

### 4. All cases together

{: align="center"}
|![All]({{ fig_loc }}ThreeScenariosDecline.png)|

Source: Authors’ calculation.<br>
Note:	See Fujii and Nakata (2021) for a detailed discussion of the scenario.
-->
