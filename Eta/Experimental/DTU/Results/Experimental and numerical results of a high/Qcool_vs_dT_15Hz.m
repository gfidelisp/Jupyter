
dT10=[2.35 13.71 18.04 18.86];
Qc10=[360.7	240.2	120.1	0];
dT15=[5.33 10.11 18.52 20.36 22.04];
Qc15=[321	270.94	180.08	90.7	0];
dT20=[5.45 16.06 21.02 24.48];
Qc20=[300.56	201.17	100.94	0];
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
plot(dT10,Qc10,'-ko','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',8)
plot(dT15,Qc15,'-bs','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',8)
plot(dT20,Qc20,'-rv','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',8)


% Create xlabel
xlabel('Temperature span (K)','FontSize',16);

% Create ylabel
ylabel('Cooling power (W)','FontSize',16);
h = legend('1 Hz','1.5 Hz','2 Hz')
set(h,'Fontsize',14)

print('-depsc','fig_coolingpower_freq.eps');
