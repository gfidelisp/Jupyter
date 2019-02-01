close all
clear variables


figure4 = figure('PaperPositionMode', 'auto');
fig4 = axes('Parent',figure4,'Layer','top','FontSize',16);
grid on
hold all
box on

data = load('B(z)_v=0til360_trin22.5.txt');

data = data(:,[1:5,7:38,40:length(data(1,:))]); %Excluding colomns 5 and 31

data(isnan(data)) = 0;

n_data_points = 3:length(data(:,1));
z = data(n_data_points,1)*10;
angles = unique(data(1,:));
angles = angles(1:length(angles)-1);

% Color_arr = [0.8, 0, 0.8,  0, 0, 1,  0, 0.5, 0,  1, 0, 0];
Color_arr = get(gcf,'DefaultAxesColorOrder');
Color_arr(10,:) = Color_arr(3,:);
Color_arr(3,:) = Color_arr(2,:);
Color_arr(2,:) = Color_arr(10,:);
Line_arr = {'-' '-.' '--' ':'};
Marker_arr = {'s' 'x' 'd'};


%--- Fake legend for fig4
for i = [1 3]
    plot(fig4,0,'k','Visible','off','Marker',Marker_arr{i},'Markersize',10,'Linestyle','none')
end
plot(fig4,0,'-','Visible','off','Color','k','Linewidth',1,'Marker','none')
plot(fig4,0,'--','Visible','off','Color','k','Linewidth',2,'Marker','none')

legend_fig4 = legend(fig4,'High field','Low field','Simulation data','Ends of magnet','Location','North');

for i = 1:length(angles)
    indx_angle = find(data(1,:) == angles(i));   
    for j = 1:3
        indx_component = find(data(2,indx_angle) == j);
        data_arr(i,j,:) = sum(data(n_data_points,indx_angle(indx_component)),2)/length(indx_component);
    end
end

for i = [1 3]
    plot(fig4,z,sqrt(squeeze(data_arr(i,1,:)).^2+squeeze(data_arr(i,2,:)).^2+squeeze(data_arr(i,3,:)).^2),'k','Marker',Marker_arr{i},'Markersize',10,'Linestyle','none')
end

plot(fig4,[-125 -125],[-1 2],'k--','Linewidth',2)
plot(fig4,[125 125],[-1 2],'k--','Linewidth',2)

norm_high   = load('norm_high.txt');
norm_low    = load('norm_low.txt');

plot(fig4,norm_high(:,1)*1000,norm_high(:,2),'k-');
plot(fig4,norm_low(:,1)*1000,norm_low(:,2),'k-');

xlim([-175 175]) 
ylim([-0.1 1.3]) 
xlabel('z [mm]')
ylabel('||B|| [T]');
legend_fig4_pos = get(legend_fig4,'Position');
set(legend_fig4,'Position',[legend_fig4_pos(1) legend_fig4_pos(2)-0.2 legend_fig4_pos(3) legend_fig4_pos(4)]);
print('-depsc','Norm_flux_z.eps');