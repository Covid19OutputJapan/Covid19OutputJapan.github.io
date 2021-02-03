---
title: "Covid-19 and Output in Japan"
keywords: sample homepage
tags: [nationwide]
sidebar: home_sidebar
permalink: index.html
summary:
---

{% assign fig_loc = "./archives/20210202/Figures/" %}
<!-- csv files must be in the "_data" folder -->
{% assign datafile_oneweek = site.data.oneweek20210126 %}
{% assign datafile_fourweek = site.data.fourweek20210126 %}

## Last update on February 3, 2021

{% include link_to_code.html %}

### 1. Conditional Projections of Covid-19

{: align="center"}
|![Projection]({{ fig_loc }}VariablesProjection.png)|

Source: Authors' calculation.<br>
Note: Red line:average output loss=1.35%. Black line:average output loss=1.65%. Blue line:average output loss=2.5%. Weekly frequency.

### 2. Projected relationship between Covid-19 and output

{: align="center"}
|![TradeoffUB]({{ fig_loc }}BaselineTradeoffUBp.png)|

Source: Authors' calculation.<br>
Note: Black line:the most recent week. Red line:one week earlier. Blue line:two weeks earlier. The vertical axis shows the number of cumulative deaths by the end of the next 12 months. The horizontal axis shows the average output loss over the next twelve months. The horizontal dashed line indicates the total number of Covid-19 deaths during 2020. The darkest and the second darkest grey areas indicate 20- and 40-percent confidence sets, respectively. The second lightest and the lightest grey areas indicated 60- and 80-percent confidence sets, respectively.

### 3. Forecast Errors

#### i. One-week horizon

<!--
{: align="center"}
<table>
  {% for row in datafile_oneweek %}
    {% if forloop.first %}
      <tr><th></th>
      <th> "Conditional"<br>forecast<br>from last week </th>
      <th> <br><br>Actual </th>
      <th> <br>"Conditional"<br>forecast error </th>
      </tr>
    {% endif %}
    <tr>
      {% for pair in row %}
        <td>
        {% if forloop.first %}
          <b>{{ pair[1] }}</b>
        {% else %}
          {% assign temp = pair[1] | plus:0 %}
          {% if temp > 0 %}
            <span style="color: black; ">{{ pair[1] }}</span>
          {% else %}
            <span style="color: red; ">{{ pair[1] }}</span>
          {% endif %}
        {% endif %}
        </td>
      {% endfor %}
    </tr>
  {% endfor %}
</table>
-->

{: align="center"}
|    | "Conditional"<br>forecast<br>from last week | <br><br>Actual | <br>"Conditional"<br>forecast error |
| ---- | ---- | ---- | ---- |
| **New Cases** | 34,281   |  24,238  | <span style="color: black; ">10,043</span> |
| **New Deaths** |   449  | 637  | <span style="color: red; ">-190</span> |

Source: Authors' calculation.<br>
Note: Our model provides projections of Covid-19 conditional on paths of output but does not provide the projection of output. Thus, we compare the actual Covid-19 outcomes with the projection of them if the model had known the realized path of output.  

<!--
#### ii. Four-week horizon

{: align="center"}
<table>
  {% for row in datafile_fourweek %}
    {% if forloop.first %}
    <tr><th></th>
    <th> "Conditional"<br>forecast<br>from 4 weeks ago </th>
    <th> <br><br>Actual </th>
    <th> <br>"Conditional"<br>forecast error </th>
    </tr>
    {% endif %}
    <tr>
      {% for pair in row %}
        <td>
        {% if forloop.first %}
          <b>{{ pair[1] }}</b>
        {% else %}
          {% assign temp = pair[1] | plus:0 %}
          {% if temp > 0 %}
            <span style="color: black; ">{{ pair[1] }}</span>
          {% else %}
            <span style="color: red; ">{{ pair[1] }}</span>
          {% endif %}
        {% endif %}
        </td>
      {% endfor %}
    </tr>
  {% endfor %}
</table>

{: align="center"}
|    | "Conditional"<br>forecast<br>from 4 weeks ago | <br><br>Actual | <br>"Conditional"<br>forecast error |
| ---- | ---- | ---- | ---- |
| **New Cases** |  83,138  |  129,454  | <span style="color: red; ">-46,315</span> |
| **New Deaths** |   1,004  |    1,459 | <span style="color: red; ">-454</span> |

Source: Authors' calculation.<br>
Note: Our model provides projections of Covid-19 conditional on paths of output but does not provide the projection of output. Thus, we compare the actual Covid-19 outcomes with the projection of them if the model had known the realized path of output.  
-->

### 4. Real-time Evaluation of the Model's Forecasting Performance

#### i. One-week horizon

{: align="center"}
|![ForecaseErrors]({{ fig_loc }}ForecastErrors.png)|

Source: Authors' calculation.<br>
Note: The red lines---"Forecast"---are what the model would have predicted given the data available up to that point.

<!--
#### New Cases

{: align="center"}
|![ForecastErrorsN]({{ fig_loc }}ForecastErrorsN.png)|

Source: Authors' calculation.<br>
Note: The red lines---"Forecast"---are what the model would have predicted given the data available up to that point.

#### New Deaths

{: align="center"}
|![ForecastErrorsD]({{ fig_loc }}ForecastErrorsD.png)|

Source: Authors' calculation.<br>
Note: The red lines---"Forecast"---are what the model would have predicted given the data available up to that point.
-->

<!--
### 5. Chart of the Week

{: align="center"}
|![TradeoffUB]({{ fig_loc }}ChartOfTheWeek.png)|

Source: Authors' calculation.<br>
Note: The solid black line and grey fan chart are the same as in Chart 2. The red line is the relationship between the output loss and the number of and cumulative deaths if the pace of vaccine distribution is twice as fast as in the baseline case.
-->
