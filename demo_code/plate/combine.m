clear
clc
force_address1='gamma_45_forces.csv';
position_address1='gamma_45_positions.csv';
force_address2='gamma_zero_forces.csv';
position_address2='gamma_zero_positions.csv';
slope_list1=get_slope(force_address1,position_address1);
slope_list2=get_slope(force_address2,position_address2);
s_x=[slope_list1(3,:);slope_list2(3,:)];
s_z=[slope_list1(4,:);slope_list2(4,:)];
figure
subplot(1,2,1)
xvalues=categorical(slope_list1(1,:));
yvalues=categorical([slope_list1(2),slope_list2(2)]);
hx=heatmap(xvalues,yvalues,s_z,'Colormap',jet);
hx.Title = '\alpha_z (N/cm^3)';
hx.XLabel = '\gamma';
hx.YLabel = '\beta';
hx.FontSize=20;
subplot(1,2,2)
hz=heatmap(xvalues,yvalues,s_x,'Colormap',jet);
hz.Title = '\alpha_x (N/cm^3)';
hz.XLabel = '\gamma';
hz.YLabel = '\beta';
hz.FontSize=20;