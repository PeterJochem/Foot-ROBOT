clc
clear
init_env();
force_address='output_plate_forces.csv';
position_address='output_plate_positions.csv';
total_data_force=csvread(force_address);
total_data_position=csvread(position_address);
slope_list=[];
for i= 1:1
    range_start=2+802*(i-1);
    range_terminal=7500+802*(i-1);
    extract_start=1502;
    extract_terminal=7450;
    range=range_start:range_terminal;
    extract_range=extract_start:extract_terminal;
    raw_data_start=1;
    Area=3.81*2.54;
    fitted_start=1;
    fitted_terminal=1;
    gamma=total_data_force(1+802*(i-1),1);
    beta=total_data_force(1+802*(i-1),2);
    %create some arrays to store the forces and positions
    time=total_data_force(range,1);
    f_x=total_data_force(range,2);
    f_y=total_data_force(range,3);
    f_z=total_data_force(range,4);
    pos_z=total_data_position(range,4);
    %create some arrays to store the forces and positions in extract range
    pos_z_extract=pos_z(extract_range);
    f_x_extract=f_x(extract_range);
    f_y_extract=f_y(extract_range);
    f_z_extract=f_z(extract_range);
    %granular materials top
    %to find the start index for raw data
    raw_data_start=find_index(raw_data_start,0.01,f_z);
    gm_top=total_data_position(300+802*(i-1),5);
    raw_data_terminal=length(time);
    raw_data_range=raw_data_start:raw_data_terminal;
    %to find the terminal index for fitting data
    fitted_start=find_index_extract(fitted_start,7.62,gm_top-pos_z_extract);
    %to find the terminal index for fitting data
    fitted_terminal=find_index_extract(fitted_terminal,2.54,gm_top-pos_z_extract);
    fitted_data_range=fitted_start:fitted_terminal;
%     fit1=find_index_extract(fitted_start,7.62,gm_top-pos_z_extract);
%     fit2=find_index_extract(fitted_start,5.08,gm_top-pos_z_extract);
%     fit3=find_index_extract(fitted_start,2.54,gm_top-pos_z_extract);
%     fitx_array=[f_x_extract(fit1);f_x_extract(fit2);f_x_extract(fit3)];
%     fitz_array=[f_z_extract(fit1);f_z_extract(fit2);f_z_extract(fit3)];
%     dep_array=[gm_top-pos_z_extract(fit1);gm_top-pos_z_extract(fit2);gm_top-pos_z_extract(fit3)];
    
    dep=gm_top-pos_z_extract(fitted_data_range);
    s_x=polyfit(dep,f_x_extract(fitted_data_range)/Area,1);
    s_z=polyfit(dep,f_z_extract(fitted_data_range)/Area,1);
%     s_x=polyfit(dep_array,fitx_array/Area,1);
%     s_z=polyfit(dep_array,fitz_array/Area,1);
    sx=polyval(s_x,dep);
    sz=polyval(s_z,dep);
    %plot figures
    subplot(3,2,i)
    plot(gm_top-pos_z(extract_range),f_x(extract_range)/Area);
    hold on
    plot(gm_top-pos_z(extract_range),f_z(extract_range)/Area);
    hold on
    plot(dep,sx,'LineWidth',4);
    hold on
    plot(dep,sz,'LineWidth',4);
    title("$\gamma$: "+ gamma + " $\beta$: "+beta);
    xlabel('depth (cm)','FontSize',15);
    ylabel('Stress (N/cm^2)','FontSize',15);
    legend('Stress x','Strees z','Fitted-stress x','Fitted-stress z','FontSize',10);
    set(gca,'FontSize',15)
    slope=[gamma;beta;s_x(1);s_z(1)];
    slope_list=[slope,slope_list];
end
hold off;
figure
subplot(1,2,1)
xvalues=categorical(slope_list(1,:));
yvalues=categorical(slope_list(2));
hx=heatmap(xvalues,yvalues,slope_list(3,:),'Colormap',jet);
hx.Title = '\alpha_x (N/cm^3)';
hx.XLabel = '\gamma';
hx.YLabel = '\beta';
caxis([-0.11 0.11])
hx.FontSize=20;
subplot(1,2,2)
hz=heatmap(xvalues,yvalues,slope_list(4,:),'Colormap',jet);
hz.Title = '\alpha_z (N/cm^3)';
hz.XLabel = '\gamma';
hz.YLabel = '\beta';
caxis([-0.25 0.25])
hz.FontSize=20;