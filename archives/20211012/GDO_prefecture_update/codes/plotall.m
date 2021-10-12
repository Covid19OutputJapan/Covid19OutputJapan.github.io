
for var_index = 1:size(Pref_Vars,2)
    plot(normalize(Data_struct.(Pref_Vars{var_index})))
    hold on
end

text_legend_1 = 'legend(';
text_legend_2 = '';
text_legend_3 = ');';
        
for var_index = 1:size(Pref_Vars,2)
    if var_index == size(Pref_Vars,2)
        text_legend_2 = strcat(text_legend_2,'''',Pref_Vars{var_index},'''')
    else
        text_legend_2 = strcat(text_legend_2,'''',Pref_Vars{var_index},'''',',')
    end
end
        
strcat(text_legend_1,text_legend_2,text_legend_3)
eval(strcat(text_legend_1,text_legend_2,text_legend_3))