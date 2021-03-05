
clear variables
close all
iPC=0;
if iPC==1
    home = '\Users\shcor\Dropbox\fujii_nakata\Website\Codes\';
else
    home = '/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/Codes/';
end
 variant = importdata('/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Data/variant/S_N501.csv');  % Import corss-country Covid-19 variant data
 Data = variant.data;
 Country = variant.textdata(1,2:end);
 Date = variant.textdata(2:end,1);
 W = datetime(Date,'InputFormat','yyyy/MM/dd');
 