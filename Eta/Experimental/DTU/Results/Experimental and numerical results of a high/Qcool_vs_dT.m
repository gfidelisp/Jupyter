
dT200=[2.29 9.69 15.11 15.73 16.39];
Qc200=[188.5	150.1	100.8	50.1	0];
dT300=[3.36	13.38 16.68 17.68];
Qc300=[300.6	200.8	101.8	0];
dT400=[2.35	13.71 18.04 18.86];
Qc400=[360.7	240.2	120.1	0];
% Qnorm=667.75751;
% Qc400=Qc400/Qnorm;
% Qc600=Qc600/Qnorm;
% Qc800=Qc800/Qnorm;
%plot results

% Create figure
figure1 = figure('PaperSize',[20.98 29.68]);

% Create axes
axes('Parent',figure1,'FontSize',14);
% Uncomment the following line to preserve the Y-limits of the axes
%ylim([0 1.05]);
box('on');
grid('on');
hold('all');

% Create plot
plot(dT200,Qc200,'-ko','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',8)
plot(dT300,Qc300,'-bs','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',8)
plot(dT400,Qc400,'-rv','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',8)


% Create xlabel
xlabel('Temperature span (K)','FontSize',16);

% Create ylabel
ylabel('Cooling power (W)','FontSize',16);
h = legend('200 L/h','300 L/h','400 L/h')
set(h,'Fontsize',14)

print('-depsc','fig_coolingpower_load.eps');
