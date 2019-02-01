close all 
clear

%% Temperature span vs. frequency

fexp = [0.25	0.5	1	1.5	2.25	3	4	5	6	7	8	9 10];
dTexp = [2.42	9.84	15.78	17.52	18.9	18.74	18.34	17.69	16.76	15.93	14.97	13.5 12.09];
fexp_400 = [0.25 0.5 1 1.5 2.25 3 4];
dTexp_400 = [3.92 11.69 15.4 16.81 17.4 16.8 15.77];
fexp_200 = [0.25 0.5 1 1.5];
dTexp_200 = [4.28 7.83 8.59 4.99];
fnoloss = [0.25	0.5	1 2 4 6 8 10];
dTnoloss = [1.02 15.51 20.64 22.01 21.89 20.83 19.28 17.31];
floss = [0.25 0.5 1 2 4 6 8 10];
dTloss = [1.04 14.53 19.16 19.96 18.91 16.88 12.05 4.61];

% Qnorm=667.75751;
% Qc400=Qc400/Qnorm;
% Qc600=Qc600/Qnorm;
% Qc800=Qc800/Qnorm;
% plot results

figure1 = figure('PaperSize',[20.98 29.68]);
axes('Parent',figure1,'FontSize',14);
%ylim([0 1.05]);
box('on');
grid('on');
hold('all');
%plot(fexp,dTexp,':ko','Linewidth',2','MarkerFaceColor','k','MarkerSize',8)
%plot(fexp_400,dTexp_400,':b^','Linewidth',2','MarkerFaceColor','b','MarkerSize',8)
%plot(fexp_200,dTexp_200,':rs','Linewidth',2','MarkerFaceColor','r','MarkerSize',8)
plot(fnoloss,dTnoloss,'--ko','Linewidth',1.5','MarkerSize',6)
plot(floss,dTloss,':kx','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',8)
ylabel('Temperature span (K)','FontSize',16);
xlabel('Frequency (Hz)','FontSize',16);
h = legend('500 L/h - Numerical','500 L/h - Num. losses','400 L/h','200 L/h','Numerical','Numerical with losses');
set(h,'Fontsize',14)
print('-depsc','fig_dT_vs_f.eps');

%% Power of motor with big pulley

freq_motor_296 = [3 4 5 6 7 8];
Power_motor_296 = [272 316 375 442 502 540];
freq_motor_dT = [0.25	0.5	1	1.5	2.25	3	4	5	6	7	8	9	10];
Power_motor_dT = [100	159	220	255	370	415	433	470	510	549	635	625	726]

figure1 = figure('PaperSize',[20.98 29.68]);
axes('Parent',figure1,'FontSize',14);
%xlim([2 9]);
box('on');
grid('on');
hold('all');
plot(freq_motor_dT,Power_motor_dT,':ko','Linewidth',2','MarkerFaceColor','k','MarkerSize',8)
plot(freq_motor_296,Power_motor_296,':ko','Linewidth',2','MarkerSize',8)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('Power of motor (W)','FontSize',16);
h = legend('Temperature span - 500 L/h','296 K - 100 L/h');
set(h,'Fontsize',14)
print('-depsc','fig_Wmotor_vs_f.eps');


%% Temperature as a function of Cooling capacity for different frequencies

Q_2 = [631	400	200	100]; 
dT_2 = [0.19	11.23	15.86	18.52];
Q_4 = [591	400	200	100];
dT_4 = [1.07	10.54	16.05	19.15];
Q_6 = [520	400	200	100];
dT_6 = [0.26	8.52	14.27	17.31];
Q_8 = [400	200	100];
dT_8 = [1.52	12.63	15.31]; 
Q_10 = [370	200	100];
dT_10 = [0.9	11.43	14.46];

Q_2n = [527.30 497.56 455.23 400.80 283.16 155.97 38.20]; 
dT_2n = [1 4 8 12 16 20 24];
Q_4n = [466.49 437.24 397.22 352.39 245.52 129.63 23.50];
dT_4n = [1 4 8 12 16 20 24];
Q_6n = [395.10 365.04 325.02 286.53 190.15 85.10 -10.95];
dT_6n = [1 4 8 12 16 20 24];
Q_8n = [320.03 288.93 248.38 215.27 130.20 36.32 -49.86];
dT_8n = [1 4 8 12 16 20 24]; 
Q_10n = [245.56 213.50 172.28 143.60 70.11 -12.38 -89.23];
dT_10n = [1 4 8 12 16 20 24];

% Qnorm=667.75751;
% Qc400=Qc400/Qnorm;
% Qc600=Qc600/Qnorm;
% Qc800=Qc800/Qnorm;
figure1 = figure('PaperSize',[20.98 29.68]);
axes('Parent',figure1,'FontSize',14);
xlim([0 700])
ylim([0 20]);
box('on');
grid('on');
hold('all');
plot(Q_2,dT_2,':k^','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',7)
plot(Q_4,dT_4,':bo','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',7)
plot(Q_6,dT_6,':rs','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',7)
plot(Q_8,dT_8,':m^','Linewidth',1.5','MarkerSize',7)
plot(Q_10,dT_10,':go','Linewidth',1.5','MarkerSize',7)

% plot(Q_2n,dT_2n,'--r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',8)
% plot(Q_4n,dT_4n,'--go','Linewidth',1.5','MarkerFaceColor','g','MarkerSize',8)
% plot(Q_6n,dT_6n,'--cs','Linewidth',1.5','MarkerFaceColor','c','MarkerSize',8)
% plot(Q_8n,dT_8n,'--k*','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',8)
% plot(Q_10n,dT_10n,'--bx','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',8)

% plot(fnoloss,dTnoloss,':k','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',8)
% plot(floss,dTloss,'--k','Linewidth',2','MarkerFaceColor','k','MarkerSize',8)
ylabel('Temperature span (K)','FontSize',14);
xlabel('Cooling capacity (W)','FontSize',14);
h = legend('2 Hz','4 Hz','6 Hz','8 Hz','10 Hz');
set(h,'Fontsize',14)
print('-depsc','fig_dT_vs_Q.eps');



%% Comparison between experiment and numerical results

Q_2 = [631	400	200	100]; 
dT_2 = [0.19	11.23	15.86	18.52];
Q_4 = [591	400	200	100];
dT_4 = [1.07	10.54	16.05	19.15];
Q_6 = [520	400	200	100];
dT_6 = [0.26	8.52	14.27	17.31];
Q_8 = [400	200	100];
dT_8 = [1.52	12.63	15.31]; 
Q_10 = [370	200	100];
dT_10 = [0.9	11.43	14.46];

Q_2n = [527.30 497.56 455.23 400.80 283.16 155.97 38.20]; 
dT_2n = [1 4 8 12 16 20 24];
Q_4n = [466.49 437.24 397.22 352.39 245.52 129.63 23.50];
dT_4n = [1 4 8 12 16 20 24];
Q_6n = [395.10 365.04 325.02 286.53 190.15 85.10 -10.95];
dT_6n = [1 4 8 12 16 20 24];
Q_8n = [320.03 288.93 248.38 215.27 130.20 36.32 -49.86];
dT_8n = [1 4 8 12 16 20 24]; 
Q_10n = [245.56 213.50 172.28 143.60 70.11 -12.38 -89.23];
dT_10n = [1 4 8 12 16 20 24];

Q_2nnl = [547.30 523.19 490.26 445.43 337.42 219.88 111.78];
dT_2nnl = [1 4 8 12 16 20 24];
Q_4nnl = [506.49 482.87 452.26 417.01 319.78 213.54 117.07];
dT_4nnl = [1 4 8 12 16 20 24];
Q_10nnl = [345.56 319.13 287.31 268.22 204.37 131.54 64.35];
dT_10nnl = [1 4 8 12 16 20 24]; 

figure3 = figure('PaperSize',[20.98 29.68]);
axes('Parent',figure3,'FontSize',12);
xlim([0 600])
%ylim([0 20]);
box('on');
grid('on');
hold('all');

plot(Q_4,dT_4,'-ko','Linewidth',2','MarkerFaceColor','k','MarkerSize',8)
plot(Q_4nnl,dT_4nnl,'--ko','Linewidth',1.5','MarkerSize',6)
plot(Q_4n,dT_4n,':kx','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',8)

ylabel('Temperature span (K)','FontSize',14);
xlabel('Cooling capacity (W)','FontSize',14);
h = legend('Experimental','Numerical','Numerical with losses');
print('-depsc','fig_dT_vs_Q_4Hz.eps');

% Create figure
figure3 = figure('PaperSize',[20.98 29.68]);
% Create axes
axes('Parent',figure3,'FontSize',12);
% Uncomment the following line to preserve the Y-limits of the axes
xlim([0 600])
%ylim([0 20]);
box('on');
grid('on');
hold('all');


figure4 = figure('PaperSize',[20.98 29.68]);
axes('Parent',figure4,'FontSize',12);
xlim([0 400])
%ylim([0 20]);
box('on');
grid('on');
hold('all');

plot(Q_10,dT_10,'-ko','Linewidth',2','MarkerFaceColor','k','MarkerSize',8)
plot(Q_10nnl,dT_10nnl,'--ko','Linewidth',1.5','MarkerSize',6)
plot(Q_10n,dT_10n,':kx','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',8)

% plot(Q_10,dT_10,'-ko','Linewidth',2','MarkerFaceColor','k','MarkerSize',8)
% plot(Q_10nnl,dT_10nnl,':','Linewidth',1.5')
% plot(Q_10n,dT_10n,'--','Linewidth',2')

ylabel('Temperature span (K)','FontSize',14);
xlabel('Cooling capacity (W)','FontSize',14);
h = legend('Experimental','Numerical','Numerical with losses');
print('-depsc','fig_dT_vs_Q_10Hz.eps');

figure5 = figure('PaperSize',[20.98 29.68]);
axes('Parent',figure5,'FontSize',12);
xlim([0 700])
%ylim([0 20]);
box('on');
grid('on');
hold('all');

plot(Q_2,dT_2,'-ko','Linewidth',2','MarkerFaceColor','k','MarkerSize',8)
plot(Q_2nnl,dT_2nnl,'--ko','Linewidth',1.5','MarkerSize',6)
plot(Q_2n,dT_2n,':kx','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',8)

% plot(Q_2,dT_2,'-ko','Linewidth',2','MarkerFaceColor','k','MarkerSize',8)
% plot(Q_2nnl,dT_2nnl,':','Linewidth',1.5')
% plot(Q_2n,dT_2n,'--','Linewidth',2')

ylabel('Temperature span (K)','FontSize',14);
xlabel('Cooling capacity (W)','FontSize',14);
h = legend('Experimental','Numerical','Numerical with losses');
print('-depsc','fig_dT_vs_Q_2Hz.eps');


%% Maximum cooling capacity as a function of frequency at Tch = 20C

f = [2	4	6	8	10];   
Qmax = [635	613	524	427	385];
fnoloss = [0.5	1 2 4 6 8 10];
Qmaxnoloss = [566 568 555 514 463 408 354];
floss = [0.5 1 2 4 6 8 10];
Qmaxloss = [563 560 538 476 405 330 256];

% Create figure
figure6 = figure('PaperSize',[20.98 29.68]);

% Create axes
axes('Parent',figure6,'FontSize',12);
%ylim([0 1.05]);
box('on');
grid('on');
hold('all');

% Create plot
plot(f,Qmax,'-ko','Linewidth',2','MarkerFaceColor','k','MarkerSize',8)
plot(fnoloss,Qmaxnoloss,'--ko','Linewidth',1.5','MarkerSize',6)
plot(floss,Qmaxloss,':kx','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',8)
% Create xlabel
ylabel('Cooling capacity (W)','FontSize',14);
% Create ylabel
xlabel('Operating frequency (Hz)','FontSize',14);
h = legend('Experimental','Numerical','Numerical with losses');
print('-depsc','fig_Qmax_vs_f.eps');