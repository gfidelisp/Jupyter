
freq=[0.25	0.5	1	1.5	2.25	3	4	5	6	7	8];
Qc=[2.42	9.84	15.78	17.52	18.9	18.74	18.34	17.69	16.76	15.93	14.97];
% Qnorm=667.75751;
% Qc400=Qc400/Qnorm;
% Qc600=Qc600/Qnorm;
% Qc800=Qc800/Qnorm;
%plot results

% Create figure
figure1 = figure('PaperSize',[20.98 29.68]);

% Create axes
axes('Parent',figure1,'FontSize',12);
% Uncomment the following line to preserve the Y-limits of the axes
%ylim([0 1.05]);
box('on');
grid('on');
hold('all');

% Create plot
plot(freq,Qc,'-ko','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',8)


% Create xlabel
ylabel('Temperature span (K)','FontSize',14);

% Create ylabel
xlabel('Operating frequency (Hz)','FontSize',14);

print('-depsc','fig_dT_highfreq.eps');
