% for csv files
clear all;
date1 = '0121_2021';
date2 = '20210121';
load(strcat('Japan_',date1,'.mat'));
nf = java.text.DecimalFormat; % for commas for each 3 digits

% one-week
oneweek_mat = cell(3,4);
oneweek_mat{1,1} = "Type"
oneweek_mat{1,2} = "Forecast"
oneweek_mat{1,3} = "Actual"
oneweek_mat{1,4} = "Error"
% if we have yet another forecast 
%oneweek_mat{1,5} = "Forecast2"
oneweek_mat{2,1} = "New Cases"
oneweek_mat{2,2} = string(nf.format(round(dNForecast(Tdata,1),0)));
oneweek_mat{2,3} = string(nf.format(round(dNActual(Tdata,1),0)));
oneweek_mat{2,4} = string(nf.format(round(dNForecast(Tdata,1)-dNActual(Tdata,1),0)));
oneweek_mat{3,1} = "New Deaths"
oneweek_mat{3,2} = string(nf.format(round(dDForecast(Tdata,1),0)));
oneweek_mat{3,3} = string(nf.format(round(dDActual(Tdata,1),0)));
oneweek_mat{3,4} = string(nf.format(round(dDForecast(Tdata,1)-dDActual(Tdata,1),0)));

writecell(oneweek_mat, strcat('oneweek',date2,'.csv'));

% four-week
fourweek_mat = cell(3,4);
fourweek_mat{1,1} = "Type"
fourweek_mat{1,2} = "Forecast"
fourweek_mat{1,3} = "Actual"
fourweek_mat{1,4} = "Error"
fourweek_mat{2,1} = "New Cases"
fourweek_mat{2,2} = string(nf.format(round(dNForecast(Tdata-3,2),0)));
fourweek_mat{2,3} = string(nf.format(round(dNActual(Tdata-3,2),0)));
fourweek_mat{2,4} = string(nf.format(round(dNForecast(Tdata-3,2)-dNActual(Tdata-3,2),0)));
fourweek_mat{3,1} = "New Deaths"
fourweek_mat{3,2} = string(nf.format(round(dDForecast(Tdata-3,2),0)));
fourweek_mat{3,3} = string(nf.format(round(dDActual(Tdata-3,2),0)));
fourweek_mat{3,4} = string(nf.format(round(dDForecast(Tdata-3,2)-dDActual(Tdata-3,2),0)));

writecell(fourweek_mat, strcat('fourweek',date2,'.csv'));
