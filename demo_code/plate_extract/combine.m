clear
clc
force_address1='beta_90_forces.csv';
position_address1='beta_90_positions.csv';
force_address2='beta_60_forces.csv';
position_address2='beta_60_positions.csv';
force_address3='beta_45_forces.csv';
position_address3='beta_45_positions.csv';
force_address4='beta_30_forces.csv';
position_address4='beta_30_positions.csv';
force_address5='beta_0_forces.csv';
position_address5='beta_0_positions.csv';
force_address6='beta_-45_forces.csv';
position_address6='beta_45_positions.csv';
slope_list1=get_slope(force_address1,position_address1);
slope_list2=get_slope(force_address2,position_address2);
slope_list3=get_slope(force_address3,position_address3);
slope_list4=get_slope(force_address4,position_address4);
slope_list5=get_slope(force_address5,position_address5);
slope_list6=get_slope(force_address6,position_address6);
s_x=[slope_list1(3,:);slope_list2(3,:);slope_list3(3,:);slope_list4(3,:);slope_list5(3,:);slope_list6(3,:);slope_list1(3,:)];
s_z=[slope_list1(4,:);slope_list2(4,:);slope_list3(4,:);slope_list4(4,:);slope_list5(4,:);slope_list6(4,:);slope_list1(4,:)];
figure
subplot(1,2,1)
xvalues=categorical(slope_list1(1,:));
yvalues=categorical([slope_list1(2),slope_list2(2),slope_list3(2),slope_list4(2),slope_list5(2),slope_list6(2),-90]);
hx=heatmap(xvalues,yvalues,s_z,'Colormap',jet);
hx.Title = '\alpha_z (N/cm^3)';
hx.XLabel = '\gamma';
hx.YLabel = '\beta';
hx.FontSize=20;
caxis([-0.25 0.25])
subplot(1,2,2)
hz=heatmap(xvalues,yvalues,s_x,'Colormap',jet);
hz.Title = '\alpha_x (N/cm^3)';
hz.XLabel = '\gamma';
hz.YLabel = '\beta';
caxis([-0.11 0.11])
hz.FontSize=20;