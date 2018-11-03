% Function created by Jaime A. Lozano in 30th of July 2012 to calculate
% and plot experimental and numerical results obtained from Kurt's 1D
% model. There is a thermal loss analysis and performance metrics of the
% prototype. 

%function [Qloss]=LossesPrototype(TC,TH,Tamb,freq);

close all
clear
load results_2012

%% Organize data from Experiments:

T_H_200W(:,12:13) = zeros;
T_H_400W(:,10:13) = zeros;
V_200(:,6:13) = zeros;
V_400(:,5:13) = zeros;
%For varying hot end in temperature at 200W:
freq(1,:)=T_H_200W(1,:); %Hz
V(1,:)=T_H_200W(2,:); %L/h
QC(1,:)=T_H_200W(3,:); %W
TChiller(1,:)=T_H_200W(4,:)+273; %K
Troom(1,:)=T_H_200W(5,:)+273; %K
THout(1,:)=T_H_200W(6,:)+273; %K 
TCout(1,:)=T_H_200W(7,:)+273; %K
THin(1,:)=T_H_200W(8,:)+273; %K
TCin(1,:)=T_H_200W(9,:)+273; %K
pHin(1,:)=T_H_200W(10,:); %bar
pHout(1,:)=T_H_200W(11,:); %bar
pC(1,:)=T_H_200W(12,:); %bar
Wmotor(1,:)=T_H_200W(13,:); %W
%For varying hot end in temperature at 400 W:
freq(2,:)=T_H_400W(1,:); %Hz
V(2,:)=T_H_400W(2,:); %L/h
QC(2,:)=T_H_400W(3,:); %W
TChiller(2,:)=T_H_400W(4,:)+273; %K
Troom(2,:)=T_H_400W(5,:)+273; %K
THout(2,:)=T_H_400W(6,:)+273; %K 
TCout(2,:)=T_H_400W(7,:)+273; %K
THin(2,:)=T_H_400W(8,:)+273; %K
TCin(2,:)=T_H_400W(9,:)+273; %K
pHin(2,:)=T_H_400W(10,:); %bar
pHout(2,:)=T_H_400W(11,:); %bar
pC(2,:)=T_H_400W(12,:); %bar
Wmotor(2,:)=T_H_400W(13,:); %W
%For varying volume flow rate at 200 W:
freq(3,:)=V_200(1,:); %Hz
V(3,:)=V_200(2,:); %L/h
QC(3,:)=V_200(3,:); %W
TChiller(3,:)=V_200(4,:)+273; %K
Troom(3,:)=V_200(5,:)+273; %K
THout(3,:)=V_200(6,:)+273; %K 
TCout(3,:)=V_200(7,:)+273; %K
THin(3,:)=V_200(8,:)+273; %K
TCin(3,:)=V_200(9,:)+273; %K
pHin(3,:)=V_200(10,:); %bar
pHout(3,:)=V_200(11,:); %bar
pC(3,:)=V_200(12,:); %bar
Wmotor(3,:)=V_200(13,:); %W
%For varying volume flow rate at 400 W:
freq(4,:)=V_400(1,:); %Hz
V(4,:)=V_400(2,:); %L/h
QC(4,:)=V_400(3,:); %W
TChiller(4,:)=V_400(4,:)+273; %K
Troom(4,:)=V_400(5,:)+273; %K
THout(4,:)=V_400(6,:)+273; %K 
TCout(4,:)=V_400(7,:)+273; %K
THin(4,:)=V_400(8,:)+273; %K
TCin(4,:)=V_400(9,:)+273; %K
pHin(4,:)=V_400(10,:); %bar
pHout(4,:)=V_400(11,:); %bar
pC(4,:)=V_400(12,:); %bar
Wmotor(4,:)=V_400(13,:); %W
%For varying Frequency at 200 W:
freq(5,:)=freq_200W(1,:); %Hz
V(5,:)=freq_200W(2,:); %L/h
QC(5,:)=freq_200W(3,:); %W
TChiller(5,:)=freq_200W(4,:)+273; %K
Troom(5,:)=freq_200W(5,:)+273; %K
THout(5,:)=freq_200W(6,:)+273; %K 
TCout(5,:)=freq_200W(7,:)+273; %K
THin(5,:)=freq_200W(8,:)+273; %K
TCin(5,:)=freq_200W(9,:)+273; %K
pHin(5,:)=freq_200W(10,:); %bar
pHout(5,:)=freq_200W(11,:); %bar
pC(5,:)=freq_200W(12,:); %bar
Wmotor(5,:)=freq_200W(13,:); %W

for i=1:length(TCout)
    for j=1:4
        if TCout(j,i) == 273
            TCout(j,i) = NaN;
            THin(j,i) = NaN;
        else
            TCout(j,i) = TCout(j,i);            
            THin(j,i) = THin(j,i);
        end
    end
end

%% Organize data from Numerical Simulation for Hot temperature variation

TH_N = [286.7 288.7 290.7 292.7 294.8 296 296.8 297.8 298.9 300.8 302.9 304.9];
dT_N = [0.1 0.5 1 2 4 6 8 10 12 14 16 18 20 22 24 26];
%files = {'TH_06mm_440L-h/PrototypeJAGdProto_59.txt'};
%for i=1:length(files)
%data = importdata(files(1),'\t',11);
TC_N = zeros(length(dT_N),length(TH_N));
for k=1:length(TH_N)
    for l=1:length(dT_N)
        TC_N(l,k) = TH_N(1,k)-dT_N(1,l);
    end
end
data = importdata('Simulations/TH_06mm_440L-h/PrototypeJAGdProto_59.txt', ' ', 10);
QC_N = zeros(length(dT_N(1,:)),length(TH_N(1,:)));
COP_N_N = zeros(length(dT_N(1,:)),length(TH_N(1,:)));
for j=1:length(data.data(:,1))
    indTH = find((data.data(j,2) == TH_N));
    %inddT = find(data.data(j,1) == (TH_N(indTH)-dT_N(:)));   %It does not need to create QC matrix
    indTC = find(data.data(j,1) == round(1000.*TC_N(:,indTH))./1000);
    QC_N(indTC,indTH) = data.data(j,5);
    COP_N_N(indTC,indTH) = data.data(j,9);
end
dTint = zeros(2,length(TH_N));
COPint_N = zeros(2,length(TH_N));
for i=1:length(TH_N)
    in = QC_N(:,i) ~= 0;
    dTint(1,i) = interp1(QC_N(in,i),dT_N(in),200);
    dTint(2,i) = interp1(QC_N(in,i),dT_N(in),400);
    COPint_N(1,i) = interp1(QC_N(in,i),COP_N_N(in,i),200);
    COPint_N(2,i) = interp1(QC_N(in,i),COP_N_N(in,i),400);
end
V_N(1,:) = data.data(:,3); %L/h
freq_N(1,:)= data.data(:,4); %Hz
% U_N(:)=T_H_06_440(:,6);
% Ploss_N(:)=T_H_06_440(:,7); %bar
%COP_N_N(1,:)=data.data(:,9);
Troom_N = 295;


% After variation in the mass flow rate profile at the 1D AMR model (change
% from 15 to 8.9 degrees) for 400 L/h:

TH_N_2 = [286.7 288.7 290.7 292.7 294.8 296 296.8 297.8 298.9 300.8 302.9 304.9];
dT_N_2 = [1 2 4 6 8 10 12 14 16 18 20 22 24 26 28];
%files = {'TH_06mm_440L-h/PrototypeJAGdProto_59.txt'};
%for i=1:length(files)
%data = importdata(files(1),'\t',11);
TC_N_2 = zeros(length(dT_N_2),length(TH_N_2));
for k=1:length(TH_N_2)
    for l=1:length(dT_N_2)
        TC_N_2(l,k) = TH_N_2(1,k)-dT_N_2(1,l);
    end
end
data_N_2 = importdata('Simulations/TH_06mm_400L-h_det/PrototypeJAGdProto_69.txt', ' ', 10);
QC_N_2 = zeros(length(dT_N_2(1,:)),length(TH_N_2(1,:)));
COP_N_2 = zeros(length(dT_N_2(1,:)),length(TH_N_2(1,:)));
for j=1:length(data_N_2.data(:,1))
    indTH_2 = find((data_N_2.data(j,2) == TH_N_2));
    %inddT_2 = find(data.data(j,1) == (TH_N(indTH)-dT_N(:)));   %It does not need to create QC matrix
    indTC_2 = find(data_N_2.data(j,1) == round(1000.*TC_N_2(:,indTH_2))./1000);
    QC_N_2(indTC_2,indTH_2) = data_N_2.data(j,5);
    COP_N_2(indTC_2,indTH_2) = data_N_2.data(j,9);
end
dTint_N_2 = zeros(2,length(TH_N_2));
COPint_N_2 = zeros(2,length(TH_N_2));
for i=1:length(TH_N_2)
    in_2 = QC_N_2(:,i) ~= 0;
    dTint_N_2(1,i) = interp1(QC_N_2(in_2,i),dT_N_2(in_2),200);
    dTint_N_2(2,i) = interp1(QC_N_2(in_2,i),dT_N_2(in_2),400);
    COPint_N_2(1,i) = interp1(QC_N_2(in_2,i),COP_N_2(in_2,i),200);
    COPint_N_2(2,i) = interp1(QC_N_2(in_2,i),COP_N_2(in_2,i),400);
end
V_N_2(1,:) = data_N_2.data(:,3); %L/h
freq_N_2(1,:)= data_N_2.data(:,4); %Hz
% U_N(:)=T_H_06_440(:,6);
% Ploss_N(:)=T_H_06_440(:,7); %bar
%COP_N_N(1,:)=data.data(:,9);
Troom_N_2 = 295;


% After variation in the mass flow rate profile at the 1D AMR model (change
% from 15 to 8.9 degrees) for 440 L/h:

TH_N_3 = [286.7 288.7 290.7 292.7 294.8 296 296.8 297.8 298.9 300.8 302.9 304.9];
dT_N_3 = [1 2 4 6 8 10 12 14 16 18 20 22 24 26 28];
%files = {'TH_06mm_440L-h/PrototypeJAGdProto_59.txt'};
%for i=1:length(files)
%data = importdata(files(1),'\t',11);
TC_N_3 = zeros(length(dT_N_3),length(TH_N_3));
for k=1:length(TH_N_3)
    for l=1:length(dT_N_3)
        TC_N_3(l,k) = TH_N_3(1,k)-dT_N_3(1,l);
    end
end
data_N_3 = importdata('Simulations/TH_06mm_440L-h_det/PrototypeJAGdProto_70.txt', ' ', 10);
QC_N_3 = zeros(length(dT_N_3(1,:)),length(TH_N_3(1,:)));
COP_N_3 = zeros(length(dT_N_3(1,:)),length(TH_N_3(1,:)));
for j=1:length(data_N_3.data(:,1))
    indTH_3 = find((data_N_3.data(j,2) == TH_N_3));
    %inddT_3 = find(data.data(j,1) == (TH_N(indTH)-dT_N(:)));   %It does not need to create QC matrix
    indTC_3 = find(data_N_3.data(j,1) == round(1000.*TC_N_3(:,indTH_3))./1000);
    QC_N_3(indTC_3,indTH_3) = data_N_3.data(j,5);
    COP_N_3(indTC_3,indTH_3) = data_N_3.data(j,9);
end
dTint_N_3 = zeros(2,length(TH_N_3));
COPint_N_3 = zeros(2,length(TH_N_3));
for i=1:length(TH_N_3)
    in_3 = QC_N_3(:,i) ~= 0;
    dTint_N_3(1,i) = interp1(QC_N_3(in_3,i),dT_N_3(in_3),200);
    dTint_N_3(2,i) = interp1(QC_N_3(in_3,i),dT_N_3(in_3),400);
    COPint_N_3(1,i) = interp1(QC_N_3(in_3,i),COP_N_3(in_3,i),200);
    COPint_N_3(2,i) = interp1(QC_N_3(in_3,i),COP_N_3(in_3,i),400);
end
V_N_3(1,:) = data_N_3.data(:,3); %L/h
freq_N_3(1,:)= data_N_3.data(:,4); %Hz
% U_N(:)=T_H_06_440(:,6);
% Ploss_N(:)=T_H_06_440(:,7); %bar
%COP_N_N(1,:)=data.data(:,9);
Troom_N_3 = 295;

%% Organize data from Numerical Simulation for Volumetric flow rate variation:
%Old simulations:
TH_V = 297.7; %K
V_V = [100 200 300 350 400 500 600 700 800];
dT_V = [0.1 1 4 8 12 14 16 18 20 22 24];
TC_V = zeros(1,length(dT_V));
for k=1:length(dT_V)
    TC_V(1,k) = TH_V-dT_V(1,k);
end
dataV = importdata('Simulations/Qdot/PrototypeJAGdProto_56.txt', ' ', 10);
QC_V = zeros(length(TC_V(1,:)),length(V_V(1,:)));
COP_N_V = zeros(length(TC_V(1,:)),length(V_V(1,:)));
for j=1:length(dataV.data(:,1))
    indV = find((dataV.data(j,3) == V_V));
    %inddT = find(dataV.data(j,1) == (TH_N(indTH)-dT_N(:)));   %It does not need to create QC matrix
    indTC_V = find(dataV.data(j,1) == round(1000.*TC_V)./1000);
    %indTC_V = find(dataV.data(j,1) == TC_V);
    QC_V(indTC_V,indV) = dataV.data(j,5);
    COP_N_V(indTC_V,indV) = dataV.data(j,9);
end
dTint_V = zeros(2,length(V_V));
COPint_V = zeros(2,length(V_V));
for i=1:length(V_V)
    in = QC_V(:,i) ~= 0;
    dTint_V(1,i) = interp1(QC_V(in,i),dT_V(in),200);
    dTint_V(2,i) = interp1(QC_V(in,i),dT_V(in),400);
    COPint_V(1,i) = interp1(QC_V(in,i),COP_N_V(in,i),200);
    COPint_V(2,i) = interp1(QC_V(in,i),COP_N_V(in,i),400);
end
freq_V(1,:)= dataV.data(:,4); %Hz
% U_N(:)=T_H_06_440(:,6);
% Ploss_N(:)=T_H_06_440(:,7); %bar
%COP_N_V(1,:)=dataV.data(:,9);

%After changing the mass flow rate profile (from 15 degrees to 8.9
%degrees):
TH_V_2 = 297.7; %K
V_V_2 = [100 150 200 250 300 350 400 450 500 600 700 800];
%V_V_2 = [200 300 350 400 500 600 700 800]; %Missing some flows in the simulation
dT_V_2 = [1 2 4 6 8 10 12 14 16 18 20 22 24 26 28];
TC_V_2 = zeros(1,length(dT_V_2));
for k=1:length(dT_V_2)
    TC_V_2(1,k) = TH_V_2-dT_V_2(1,k);
end
dataV_2 = importdata('Simulations/Qdot_06mm_det/Qdot_06mm.txt', ' ', 10);
QC_V_2 = zeros(length(TC_V_2(1,:)),length(V_V_2(1,:)));
COP_N_V_2 = zeros(length(TC_V_2(1,:)),length(V_V_2(1,:)));
for j=1:length(dataV_2.data(:,1))
    indV_2 = find((dataV_2.data(j,3) == V_V_2));
    %inddT = find(dataV.data(j,1) == (TH_N(indTH)-dT_N(:)));   %It does not need to create QC matrix
    indTC_V_2 = find(dataV_2.data(j,1) == round(1000.*TC_V_2)./1000);
    %indTC_V_2 = find(dataV.data(j,1) == TC_V_2);
    QC_V_2(indTC_V_2,indV_2) = dataV_2.data(j,5);
    COP_N_V_2(indTC_V_2,indV_2) = dataV_2.data(j,9);
end
dTint_V_2 = zeros(2,length(V_V_2));
COPint_V_2 = zeros(2,length(V_V_2));
for i=1:length(V_V_2)
    in = QC_V_2(:,i) ~= 0;
    dTint_V_2(1,i) = interp1(QC_V_2(in,i),dT_V_2(in),200); %spline to interpolate negative values
    dTint_V_2(2,i) = interp1(QC_V_2(in,i),dT_V_2(in),400);
    COPint_V_2(1,i) = interp1(QC_V_2(in,i),COP_N_V_2(in,i),200);
    COPint_V_2(2,i) = interp1(QC_V_2(in,i),COP_N_V_2(in,i),400);
end
%dTint_V_2(2,6) = -3.69; % extrapolate from QC_V_2 results
%dTint_V_2(2,6) = 0.37; % extrapolate from dTint_V_2 results
dTint_V_2(2,6) = (-3.69 + 0.37) / 2

freq_V_2(1,:)= dataV_2.data(:,4); %Hz
% U_N(:)=T_H_06_440(:,6);
% Ploss_N(:)=T_H_06_440(:,7); %bar
%COP_N_V(1,:)=dataV.data(:,9);

%% Organize data from Numerical Simulation for Volumetric flow rate variation at different particle sizes:
% for sphere diameter = 0.45 mm
TH_V_045 = 297.7; %K
V_V_045 = [100 150 200 250 300 350 400 450 500 600 700 800];
dT_V_045 = [1 2 4 6 8 10 12 14 16 18 20 22 24 26 28];
TC_V_045 = zeros(1,length(dT_V_045));
for k=1:length(dT_V_045)
    TC_V_045(1,k) = TH_V_045-dT_V_045(1,k);
end
dataV_045 = importdata('Simulations/Qdot_045mm_det/Qdot_045mm.txt', ' ', 10);
QC_V_045 = zeros(length(TC_V_045(1,:)),length(V_V_045(1,:)));
COP_N_V_045 = zeros(length(TC_V_045(1,:)),length(V_V_045(1,:)));
for j=1:length(dataV_045.data(:,1))
    indV_045 = find((dataV_045.data(j,3) == V_V_045));
    %inddT = find(dataV.data(j,1) == (TH_N(indTH)-dT_N(:)));   %It does not need to create QC matrix
    indTC_V_045 = find(dataV_045.data(j,1) == round(1000.*TC_V_045)./1000);
    %indTC_V_045 = find(dataV.data(j,1) == TC_V_045);
    QC_V_045(indTC_V_045,indV_045) = dataV_045.data(j,5);
    COP_N_V_045(indTC_V_045,indV_045) = dataV_045.data(j,9);
end
dTint_V_045 = zeros(2,length(V_V_045));
COPint_V_045 = zeros(2,length(V_V_045));
for i=1:length(V_V_045)
    in = QC_V_045(:,i) ~= 0;
    dTint_V_045(1,i) = interp1(QC_V_045(in,i),dT_V_045(in),200); %spline to interpolate negative values
    dTint_V_045(2,i) = interp1(QC_V_045(in,i),dT_V_045(in),400);
    COPint_V_045(1,i) = interp1(QC_V_045(in,i),COP_N_V_045(in,i),200);
    COPint_V_045(2,i) = interp1(QC_V_045(in,i),COP_N_V_045(in,i),400);
end
freq_V_045(1,:)= dataV_045.data(:,4); %Hz

% for sphere diameter = 0.5 mm
TH_V_05 = 297.7; %K
V_V_05 = [100 150 200 250 300 350 400 450 500 600 700 800];
dT_V_05 = [1 2 4 6 8 10 12 14 16 18 20 22 24 26 28];
TC_V_05 = zeros(1,length(dT_V_05));
for k=1:length(dT_V_05)
    TC_V_05(1,k) = TH_V_05-dT_V_05(1,k);
end
dataV_05 = importdata('Simulations/Qdot_05mm_det/PrototypeJAGdProto_76.txt', ' ', 10);
QC_V_05 = zeros(length(TC_V_05(1,:)),length(V_V_05(1,:)));
COP_N_V_05 = zeros(length(TC_V_05(1,:)),length(V_V_05(1,:)));
for j=1:length(dataV_05.data(:,1))
    indV_05 = find((dataV_05.data(j,3) == V_V_05));
    %inddT = find(dataV.data(j,1) == (TH_N(indTH)-dT_N(:)));   %It does not need to create QC matrix
    indTC_V_05 = find(dataV_05.data(j,1) == round(1000.*TC_V_05)./1000);
    %indTC_V_05 = find(dataV.data(j,1) == TC_V_05);
    QC_V_05(indTC_V_05,indV_05) = dataV_05.data(j,5);
    COP_N_V_05(indTC_V_05,indV_05) = dataV_05.data(j,9);
end
dTint_V_05 = zeros(2,length(V_V_05));
COPint_V_05 = zeros(2,length(V_V_05));
for i=1:length(V_V_05)
    in = QC_V_05(:,i) ~= 0;
    dTint_V_05(1,i) = interp1(QC_V_05(in,i),dT_V_05(in),200); %spline to interpolate negative values
    dTint_V_05(2,i) = interp1(QC_V_05(in,i),dT_V_05(in),400);
    COPint_V_05(1,i) = interp1(QC_V_05(in,i),COP_N_V_05(in,i),200);
    COPint_V_05(2,i) = interp1(QC_V_05(in,i),COP_N_V_05(in,i),400);
end
freq_V_05(1,:)= dataV_05.data(:,4); %Hz

% for sphere diameter = 0.7 mm
TH_V_07 = 297.7; %K
V_V_07 = [100 150 200 250 300 350 400 450 500 600 700 800];
dT_V_07 = [1 2 4 6 8 10 12 14 16 18 20 22 24 26 28];
TC_V_07 = zeros(1,length(dT_V_07));
for k=1:length(dT_V_07)
    TC_V_07(1,k) = TH_V_07-dT_V_07(1,k);
end
dataV_07 = importdata('Simulations/Qdot_07mm_det/PrototypeJAGdProto_77.txt', ' ', 10);
QC_V_07 = zeros(length(TC_V_07(1,:)),length(V_V_07(1,:)));
COP_N_V_07 = zeros(length(TC_V_07(1,:)),length(V_V_07(1,:)));
for j=1:length(dataV_07.data(:,1))
    indV_07 = find((dataV_07.data(j,3) == V_V_07));
    %inddT = find(dataV.data(j,1) == (TH_N(indTH)-dT_N(:)));   %It does not need to create QC matrix
    indTC_V_07 = find(dataV_07.data(j,1) == round(1000.*TC_V_07)./1000);
    %indTC_V_07 = find(dataV.data(j,1) == TC_V_07);
    QC_V_07(indTC_V_07,indV_07) = dataV_07.data(j,5);
    COP_N_V_07(indTC_V_07,indV_07) = dataV_07.data(j,9);
end
dTint_V_07 = zeros(2,length(V_V_07));
COPint_V_07 = zeros(2,length(V_V_07));
for i=1:length(V_V_07)
    in = QC_V_07(:,i) ~= 0;
    dTint_V_07(1,i) = interp1(QC_V_07(in,i),dT_V_07(in),200); %spline to interpolate negative values
    dTint_V_07(2,i) = interp1(QC_V_07(in,i),dT_V_07(in),400);
    COPint_V_07(1,i) = interp1(QC_V_07(in,i),COP_N_V_07(in,i),200);
    COPint_V_07(2,i) = interp1(QC_V_07(in,i),COP_N_V_07(in,i),400);
end
freq_V_07(1,:)= dataV_07.data(:,4); %Hz

% for sphere diameter = 0.75 mm
TH_V_075 = 297.7; %K
V_V_075 = [100 150 200 250 300 350 400 450 500 600 700 800];
dT_V_075 = [1 2 4 6 8 10 12 14 16 18 20 22 24 26 28];
TC_V_075 = zeros(1,length(dT_V_075));
for k=1:length(dT_V_075)
    TC_V_075(1,k) = TH_V_075-dT_V_075(1,k);
end
dataV_075 = importdata('Simulations/Qdot_075mm_det/Qdot_075mm.txt', ' ', 10);
QC_V_075 = zeros(length(TC_V_075(1,:)),length(V_V_075(1,:)));
COP_N_V_075 = zeros(length(TC_V_075(1,:)),length(V_V_075(1,:)));
for j=1:length(dataV_075.data(:,1))
    indV_075 = find((dataV_075.data(j,3) == V_V_075));
    %inddT = find(dataV.data(j,1) == (TH_N(indTH)-dT_N(:)));   %It does not need to create QC matrix
    indTC_V_075 = find(dataV_075.data(j,1) == round(1000.*TC_V_075)./1000);
    %indTC_V_075 = find(dataV.data(j,1) == TC_V_075);
    QC_V_075(indTC_V_075,indV_075) = dataV_075.data(j,5);
    COP_N_V_075(indTC_V_075,indV_075) = dataV_075.data(j,9);
end
dTint_V_075 = zeros(2,length(V_V_075));
COPint_V_075 = zeros(2,length(V_V_075));
for i=1:length(V_V_075)
    in = QC_V_075(:,i) ~= 0;
    dTint_V_075(1,i) = interp1(QC_V_075(in,i),dT_V_075(in),200); %spline to interpolate negative values
    dTint_V_075(2,i) = interp1(QC_V_075(in,i),dT_V_075(in),400);
    COPint_V_075(1,i) = interp1(QC_V_075(in,i),COP_N_V_075(in,i),200);
    COPint_V_075(2,i) = interp1(QC_V_075(in,i),COP_N_V_075(in,i),400);
end
freq_V_075(1,:)= dataV_075.data(:,4); %Hz

%% Organize data from Numerical Simulation for Frequency variation

TH_F = 298; %K
freq_F = [0.5 1 2 4 6 8 10]; %Hz
TC_F = [297 294 290 286 282 278 274];
V_F = 500;
dT_F = TH_F - TC_F;
%TC_V = zeros(1,length(dT_V));
dataF = importdata('Simulations/freq-06-298/PrototypeJAGdProto_35.txt', ' ', 10);
QC_F = zeros(length(TC_F(1,:)),length(freq_F(1,:)));
COP_N_F = zeros(length(TC_F(1,:)),length(freq_F(1,:)));
for j=1:length(dataF.data(:,1))
    indF = find((dataF.data(j,4) == freq_F));
    indTC_F = find(dataF.data(j,1) == round(1000.*TC_F)./1000);
    QC_F(indTC_F,indF) = dataF.data(j,5);
    COP_N_F(indTC_F,indF) = dataF.data(j,9);
end
dTint_F = zeros(2,length(freq_F));
COPint_F = zeros(2,length(freq_F));
for i=1:length(freq_F)
    in = QC_F(:,i) ~= 0;
    dTint_F(1,i) = interp1(QC_F(in,i),dT_F(in),200);
    dTint_F(2,i) = interp1(QC_F(in,i),dT_F(in),400);
    COPint_F(1,i) = interp1(QC_F(in,i),COP_N_F(in,i),200);
    COPint_F(2,i) = interp1(QC_F(in,i),COP_N_F(in,i),400);
end
% U_N(:)=T_H_06_440(:,6);
% Ploss_N(:)=T_H_06_440(:,7); %bar
% COP_N_F(1,:)=dataV.data(:,9);

%% Thermal losses in the TUBES and FILTER:

%Thermal resistance of the losses on the tubes and the filter at the cold end by
%conduction, free convection and radiation. Convection inside the pipes and filter is
%ignored. 

%Dimmmensions:
R_tube = 0.0125; %m
r_tube = 0.009; %m
D_cyl = 2 * R_tube; %m
L_tubes_Cu = 0.32*2; %m
L_tubes_plas = 0.4*2; %m 
A_tubes_Cu = 2*pi*R_tube*L_tubes_Cu; %m2
A_tubes_plas = 2*pi*R_tube*L_tubes_plas; %m2
t_tube = (R_tube - r_tube); %m
t_ins = 0.01; %m
R_filter = 0.045; %m
r_filter = 0.041; %m
D_filter = 2 * R_filter;
L_filter = 0.28; %m
A_filter = 2*pi*R_filter*L_filter;

%Data:
%T_s = 295;
%T_inf = 280;
%dT = T_s-T_inf;
%Properties of Air:
g = 9.8; %m/s2
beta = 3.33e-3; 
k= 0.026; %W/mK
k_Cu = 220; %W/mK
k_ins = 0.029; %W/mk
k_plas = 0.24; %W/mK
v = 15.89e-6;
alpha = 22.5e-6;
Pr = 0.707;
sigma_HT = 5.67e-8; %[W/m^2]	%radiation
epsilon = 0.85;	%radiation

% Resistances for conduction:
Rcond_tubes_Cu = t_tube / (k_Cu * A_tubes_Cu); %0.00032 K/W
Rcond_tubes_Cu_insul = t_ins / (k_ins * A_tubes_Cu); %0.00032 K/W

Rcond_tubes_plas = t_tube / (k_plas * A_tubes_plas); %0.23 K/W
Rcond_filter = (R_filter - r_filter) / (k_plas * A_filter); %0.21 K/W

% Resistances for convection: 
% It is calculated a Ra_D for each experiment (different TCout):
for i=1:length(TCout)
    for j=1:5
    if TCout(j,i) < Troom(j,i)
        Ra_D(j,i) = (g * beta .* (Troom(j,i) - TCout(j,i)) .* D_cyl^3) / (v * alpha);
        Ra_D_filter(j,i) = (g * beta * (Troom(j,i) - TCout(j,i)) * D_filter^3) / (v * alpha);
    else
        Ra_D(j,i) = (g * beta .* (TCout(j,i) - Troom(j,i)) .* D_cyl^3) / (v * alpha);
        Ra_D_filter(j,i) = (g * beta * (TCout(j,i) - Troom(j,i)) * D_filter^3) / (v * alpha);
    end
    end
end
Nus_D = (0.6 +((0.387 * Ra_D.^(1/6))/(1+(0.559/Pr)^(9/16))^(8/27))).^2;
Nus_D_filter = (0.6 +((0.387 * Ra_D_filter.^(1/6))/(1+(0.559/Pr)^(9/16))^(8/27))).^2;

h_tubes_Cu = Nus_D * k / D_cyl;  %h_tubes_Cu = 0.7154; %Kurt
h_tubes_plas = Nus_D * k / D_cyl;  %h_tubes_plas = 0.7; %Kurt
h_filter = Nus_D_filter * k / D_filter;

Rconv_tubes_Cu = 1 ./ (h_tubes_Cu * A_tubes_Cu); %~4.0 K/W
Rconv_tubes_plas = 1 ./ (h_tubes_plas * A_tubes_plas); %~3.2 K/W
Rconv_filter = 1 ./ (h_filter * A_filter); %~4.0 K/W

% Resistances for radiation:
epsilon_Cu = 0.65;
epsilon_plas = 0.91;

hrad_tubes_Cu = epsilon_Cu * sigma_HT  * (Troom.^2 + TCout.^2) .* (Troom + TCout);
hrad_tubes_plas = epsilon_plas * sigma_HT  * (Troom.^2 + TCout.^2) .* (Troom + TCout);
hrad_filter = epsilon_plas * sigma_HT  * (Troom.^2 + TCout.^2) .* (Troom + TCout);

Rrad_tubes_Cu = 1 ./ (hrad_tubes_Cu * A_tubes_Cu); %~5.5 K/W
Rrad_tubes_plas = 1 ./ (hrad_tubes_plas * A_tubes_plas); %~3.0 K/W
Rrad_filter = 1 ./ (hrad_filter * A_filter); %~2.5 K/W

%Total Resistances:
Rtubes_Cu = Rcond_tubes_Cu + Rcond_tubes_Cu_insul + ((Rconv_tubes_Cu .* Rrad_tubes_Cu)./(Rconv_tubes_Cu + Rrad_tubes_Cu));
Rtubes_plas = Rcond_tubes_plas + ((Rconv_tubes_plas .* Rrad_tubes_plas)./(Rconv_tubes_plas + Rrad_tubes_plas));
Rfilter = Rcond_filter + ((Rconv_filter .* Rrad_filter)./(Rconv_filter + Rrad_filter));

%Average Resistances:
Rtubes_Cu_av = nanmean(nanmean(Rtubes_Cu)); % = 16.07 K/W
Rtubes_plas_av = nanmean(nanmean(Rtubes_plas)); % = 1.84 K/W
Rfilter_av = nanmean(nanmean(Rfilter)); % = 1.66 K/W

%Heat loss calculations:
Qloss_tubes_Cu=(Troom-TCout)./Rtubes_Cu;  %W Heat loss in the Cu tubes
Qloss_tubes_plas=(Troom-TCout)./Rtubes_plas;  %W Heat loss in the plastic tubes
Qloss_tubes = Qloss_tubes_Cu + Qloss_tubes_plas;  %W Heat loss in the plastic tubes
Qloss_filter=(Troom-TCout)./Rfilter; %W Heat loss in the filter

%% Thermal losses in REGENERATOR BEDS: 

% It is considered to have forced convection over the regenerator beds
% since they are rotating inside the concentric magnets. It is assumed to 
% have a 'Couette flow' as the air gap is much smaller than the diameter of 
% the cylinders. So it is assumed to have a flat plate convection  where the
% Nusselt number is 7.54. The other calculated resistances are convection 
% between the fluid in the bed to the wall and thermal conduction through
% the regenerator housing wall. It is very important to notice that there 
% is a temperature gradient inside the beds. Because of its relatively high 
% thermal conductivity the temperature of the magnet and the air in between is considered to 
% have the a constant temperature equal to that of the room. For this loss
% analysis of the regenerator it is considered that the bottom and lateral
% walls of the bed are adiabatic and only heat flux occurs in the upper
% wall. 

% Dimensions:
L_bed = 0.1; %m
w_bed = 0.0125; %m
h_bed = 0.0186; %m
t_bed = 0.0016; %m Thickness of housing wall
t_air = 0.0015; %m Thickness of air layer between magnet and AMRs
A_reg = L_bed * w_bed * 24; %m
k_nylon = 0.192; %W/mK
k_air = 0.026; %W/mK 

% Thermal resistances (only in the bed length, not including the channels)
Rcond_reg = t_bed / (k_nylon * A_reg); % = 0.28 K/W
Nus_reg_out = 7.54;
h_reg_out = Nus_reg_out * k_air / t_air;
Rconv_reg_out = 1 / (h_reg_out * A_reg); % = 0.26 K/W
h_reg_in = 56.5/A_reg; %Calculated by Kurt (model)
Rconv_reg_in = 1 / (h_reg_in * A_reg); % equals to 0.018 K/W
Rrad_reg = 1 ./ (hrad_tubes_plas * A_reg); 
Rrad_reg_av = nanmean(nanmean(Rrad_reg)); %~ 6.66 K/W
%Rreg = Rcond_reg + Rconv_reg_in + Rconv_reg_out;  % = 0.55 K/W
Rreg = Rcond_reg + Rconv_reg_in + ((Rconv_reg_out * Rrad_reg_av)/(Rconv_reg_out + Rrad_reg_av));  % = 0.54 K/W

% The length of the regenerator that will be heated by ambient can be
% calculated as L_cool = ((Troom-TCout)*L_bed)/dT_bed. The heating losses across
% the regenerator are then: 
% Qloss_bed = integral(0,L_cool, (Troom-T(x))/R', dx)
% Where T(x) = dT_bed*x/L_bed, 
% Qloss_bed = (L_bed*(Troom-TCout).^2)/(2*R'*dT_bed)
% Where R'/L_bed = Rreg

% Calculate heat loss on the regenerator from experiment data: 
Qloss_reg=zeros(5,length(TCout));
for j=1:length(TCout)
    for i=1:5
        if TCout(i,j)<Troom(i,j)
            Qloss_reg(i,j) = (Troom(i,j)-TCout(i,j))^2 / (2*Rreg*(THout(i,j)-TCout(i,j)));  %from the integral
        else
            Qloss_reg(i,j) = (Troom(i,j)-TCout(i,j))/Rreg; %basically the regenerator is rejecting heat to the room
        end
    end
end

%% Calculated thermal resistances of Flowhead and Valves:

Rflowhead=1.05; %K/W, assumes non insulation on the flow head (from 2D Comsol)
Qloss_flowhead=(Troom-TCout)/Rflowhead; %W
Tvalve=10; %W/Hz --> VERIFY THIS VALUE (excel data) 
%Tvalve = 13.2 * freq - 2.5; 
%Tvalve = 2 + (10-2)*((TCout-273)./(Troom-273)); %Suppose Tvalves = 2 W/Hz at 0ºC
%Tvalve = ((2.53*freq + 6.21).^2 + (3.19*TCout -872.66).^2); %From temperature span vs. frequency experiments. Predicted - losses (without valve loss) - experiment. 
%Tvalve = ((2.53*freq + 6.21)./2 + (3.19*TCout -872.66)./2);
%Tvalve = ((2.53*freq + 6.21));
%for i=1:length(freq(:,1))
%    for j=1:length(freq(1,:))
%         indnan = isnan(Tvalve(i,j));
%         if indnan == 1
%             Tvalve(i,j) = nan;
%         else
%             %Tvalve(i,j) = Tvalve(i,j)^0.5;
%             Tvalve(i,j) = Tvalve(i,j);
%         end
%     end
% end

%Tvalve_av = nanmean(nanmean(Tvalve)); % (change with temperature?) (Freq. operation?) frequency of operation
%Qloss_valve=Tvalve.*freq; % Be careful with the frequency employed to calculate Qvalves
Qloss_valve=13.2*freq-2.5;
Qloss_tot_exp = Qloss_tubes + Qloss_flowhead + Qloss_reg + Qloss_valve + Qloss_filter;

%Calculate losses with insulation: 

%% Calculate Losses at the Hot end:

% Resistances for conduction:
Rcond_tubes_Cu_H = t_tube / (k_Cu * A_tubes_Cu/2); %0.000633 K/W It is just 1 copper tube (from THin to TH)

% Resistances for convection: 
% It is calculated a Ra_D for each experiment (different TCout):
for i=1:length(THin)
    for j=1:5
        if THin(j,i) > Troom(j,i)
            Ra_D_H(j,i) = (g * beta .* (THin(j,i) - Troom(j,i)) .* D_cyl^3) / (v * alpha);
        else
            Ra_D_H(j,i) = (g * beta .* (Troom(j,i) - THin(j,i)) .* D_cyl^3) / (v * alpha);
        end
    end
end

Nus_D_H = (0.6 +((0.387 * Ra_D_H.^(1/6))/(1+(0.559/Pr)^(9/16))^(8/27))).^2;
h_tubes_Cu_H = Nus_D_H * k / D_cyl;  %h_tubes_Cu = 0.7154; %Kurt
Rconv_tubes_Cu_H = 1 ./ (h_tubes_Cu_H * A_tubes_Cu/2); %~10.0 K/W
% Resistances for radiation:
epsilon_Cu = 0.65;
epsilon_plas = 0.91;
hrad_tubes_Cu_H = epsilon_Cu * sigma_HT  * (Troom.^2 + THin.^2) .* (Troom + THin);
Rrad_tubes_Cu_H = 1 ./ (hrad_tubes_Cu_H * A_tubes_Cu/2); %~10 K/W

%Total Resistances:
Rtubes_Cu_H = Rcond_tubes_Cu_H + ((Rconv_tubes_Cu_H .* Rrad_tubes_Cu_H)./(Rconv_tubes_Cu_H + Rrad_tubes_Cu_H));
%Average Resistances:
Rtubes_Cu_H_av = nanmean(nanmean(Rtubes_Cu_H)); % = 2.35 K/W
%Heat loss calculations for the Copper tube:
Qloss_tube_Cu_H=(THin - Troom)./Rtubes_Cu_H;  %W Heat loss in the Cu tubes
Qloss_tube_H = Qloss_tube_Cu_H;  %W 
Qloss_flowhead_H = (THin - Troom)/Rflowhead; %W
Qloss_valve_H=13.2*freq-5;
Qloss_tot_exp_H = - Qloss_tube_H - Qloss_flowhead_H  + Qloss_valve_H;
% Qloss_reg  -- What about the 'loss' in the regenerator at the hot end? 

%% Other losses

%Torque_mag = calculate from power measurements
%Eddy currents = Calculate from block and single spheres
%they depend on frequency and temperature span
% Other losses can be: Thermodynamic losses on AMR (Heat transfer)
% Viscous dissipation losses --> pump
% Axial conduction

%% Calculate losses from Numerical results for hot temperature variation

Qloss_tubes_Cu_N=zeros(length(TC_N),length(TH_N)); %K
Qloss_tubes_plas_N=zeros(length(TC_N),length(TH_N)); %K
Qloss_flowhead_N=zeros(length(TC_N),length(TH_N)); %K
Qloss_filter_N=zeros(length(TC_N),length(TH_N)); %K
Qloss_valve_N=zeros(1,length(TH_N)); %K
Qloss_reg_N=zeros(length(TC_N),length(TH_N)); %K
Qloss_tot_N=zeros(length(TC_N),length(TH_N)); %K
 
for i=1:length(TC_N)
    for j=1:length(TH_N)
        Qloss_tubes_Cu_N(i,j)=(Troom_N-TC_N(i,j))/Rtubes_Cu_av; %W
        Qloss_tubes_plas_N(i,j)=(Troom_N-TC_N(i,j))/Rtubes_plas_av; %W
        Qloss_tubes_N(i,j)=Qloss_tubes_Cu_N(i,j) + Qloss_tubes_plas_N(i,j); %W
        Qloss_flowhead_N(i,j)=(Troom_N-TC_N(i,j))/Rflowhead; %W
        Qloss_filter_N(i,j)=(Troom_N-TC_N(i,j))/Rfilter_av; %W
        
        
        %idx= searchclosest(TCout,TC_N(i,j)) %index of closest value
        %idx= nearestNeigbour(TCout,TC_N(i,j)) %index of closest value
        %ind = interp1(myCol, 1:length(myCol), [40 65 130 201], 'nearest', 'extrap')
        
        %[val,indTC] = min(abs(TCout(:)-TC_N(i,j)));
        %Qvalves_N(1,j)=Tvalves(indTC).*freq_N(:,j);
        %Qloss_valve_N(1,j)=Tvalve_av.*freq_N(:,j);  %TAKE CARE HERE! DO IT AGAIN!
        Qloss_valve_N(1,j)=13.2*freq_N(:,j)-5;
        
        if TC_N(i,j)<Troom_N
            Qloss_reg_N(i,j)=(Troom_N-TC_N(i,j))^2 / (2*Rreg*(TH_N(1,j)-TC_N(i,j)));
        else
            Qloss_reg_N(i,j)=(Troom_N-TC_N(i,j))/Rreg;
        end
    Qloss_tot_N(i,j)=Qloss_tubes_N(i,j)+Qloss_flowhead_N(i,j)+Qloss_reg_N(i,j)+Qloss_valve_N(1,j)+Qloss_filter_N(i,j);  %total heat losses
    QCloss_tot_N(i,j) = QC_N(i,j) - Qloss_tot_N(i,j);  %Cooling capacity simulated minus the losses
    end
end

% Losses for the simulations after changing mass flow rate profile
% (400L/h):
Qloss_tubes_Cu_N_2=zeros(length(TC_N_2),length(TH_N_2)); %K
Qloss_tubes_plas_N_2=zeros(length(TC_N_2),length(TH_N_2)); %K
Qloss_flowhead_N_2=zeros(length(TC_N_2),length(TH_N_2)); %K
Qloss_filter_N_2=zeros(length(TC_N_2),length(TH_N_2)); %K
Qloss_valve_N_2=zeros(1,length(TH_N_2)); %K
Qloss_reg_N_2=zeros(length(TC_N_2),length(TH_N_2)); %K
Qloss_tot_N_2=zeros(length(TC_N_2),length(TH_N_2)); %K
 
for i=1:length(TC_N_2)
    for j=1:length(TH_N_2)
        Qloss_tubes_Cu_N_2(i,j)=(Troom_N_2-TC_N_2(i,j))/Rtubes_Cu_av; %W
        Qloss_tubes_plas_N_2(i,j)=(Troom_N_2-TC_N_2(i,j))/Rtubes_plas_av; %W
        Qloss_tubes_N_2(i,j)=Qloss_tubes_Cu_N_2(i,j) + Qloss_tubes_plas_N_2(i,j); %W
        Qloss_flowhead_N_2(i,j)=(Troom_N_2-TC_N_2(i,j))/Rflowhead; %W
        Qloss_filter_N_2(i,j)=(Troom_N_2-TC_N_2(i,j))/Rfilter_av; %W
                
        %idx= searchclosest(TCout,TC_N(i,j)) %index of closest value
        %idx= nearestNeigbour(TCout,TC_N(i,j)) %index of closest value
        %ind = interp1(myCol, 1:length(myCol), [40 65 130 201], 'nearest', 'extrap')
        %[val,indTC] = min(abs(TCout(:)-TC_N(i,j)));
        %Qvalves_N(1,j)=Tvalves(indTC).*freq_N(:,j);
        %Qloss_valve_N(1,j)=Tvalve_av.*freq_N(:,j);  %TAKE CARE HERE! DO IT AGAIN!
        
        Qloss_valve_N_2(1,j)=13.2*freq_N_2(:,j)-2.5;
        
        if TC_N_2(i,j)<Troom_N_2
            Qloss_reg_N_2(i,j)=(Troom_N_2-TC_N_2(i,j))^2 / (2*Rreg*(TH_N_2(1,j)-TC_N_2(i,j)));
        else
            Qloss_reg_N_2(i,j)=(Troom_N_2-TC_N_2(i,j))/Rreg;
        end
        %Qloss_reg_N_2(1,1) = 70;
        Qloss_reg_N_2(1,1) = (294.3 - TC_N_2(1,1))^2 / (2*Rreg*(TH_N_2(1,1) - TC_N_2(1,1))); %with the actual room temperature from the experiment
    Qloss_tot_N_2(i,j)=Qloss_tubes_N_2(i,j)+Qloss_flowhead_N_2(i,j)+Qloss_reg_N_2(i,j)+Qloss_valve_N_2(1,j)+Qloss_filter_N_2(i,j);  %total heat losses
    QCloss_tot_N_2(i,j) = QC_N_2(i,j) - Qloss_tot_N_2(i,j);  %Cooling capacity simulated minus the losses
    end
end

% Losses for the simulations after changing mass flow rate profile
% (440L/h):
Qloss_tubes_Cu_N_3=zeros(length(TC_N_3),length(TH_N_3)); %K
Qloss_tubes_plas_N_3=zeros(length(TC_N_3),length(TH_N_3)); %K
Qloss_flowhead_N_3=zeros(length(TC_N_3),length(TH_N_3)); %K
Qloss_filter_N_3=zeros(length(TC_N_3),length(TH_N_3)); %K
Qloss_valve_N_3=zeros(1,length(TH_N_3)); %K
Qloss_reg_N_3=zeros(length(TC_N_3),length(TH_N_3)); %K
Qloss_tot_N_3=zeros(length(TC_N_3),length(TH_N_3)); %K
 
for i=1:length(TC_N_3)
    for j=1:length(TH_N_3)
        Qloss_tubes_Cu_N_3(i,j)=(Troom_N_3-TC_N_3(i,j))/Rtubes_Cu_av; %W
        Qloss_tubes_plas_N_3(i,j)=(Troom_N_3-TC_N_3(i,j))/Rtubes_plas_av; %W
        Qloss_tubes_N_3(i,j)=Qloss_tubes_Cu_N_3(i,j) + Qloss_tubes_plas_N_3(i,j); %W
        Qloss_flowhead_N_3(i,j)=(Troom_N_3-TC_N_3(i,j))/Rflowhead; %W
        Qloss_filter_N_3(i,j)=(Troom_N_3-TC_N_3(i,j))/Rfilter_av; %W
                
        %idx= searchclosest(TCout,TC_N(i,j)) %index of closest value
        %idx= nearestNeigbour(TCout,TC_N(i,j)) %index of closest value
        %ind = interp1(myCol, 1:length(myCol), [40 65 130 201], 'nearest', 'extrap')
        %[val,indTC] = min(abs(TCout(:)-TC_N(i,j)));
        %Qvalves_N(1,j)=Tvalves(indTC).*freq_N(:,j);
        %Qloss_valve_N(1,j)=Tvalve_av.*freq_N(:,j);  %TAKE CARE HERE! DO IT AGAIN!
        
        Qloss_valve_N_3(1,j)=13.2*freq_N_3(:,j)-5;
        
        if TC_N_3(i,j)<Troom_N_3
            Qloss_reg_N_3(i,j)=(Troom_N_3-TC_N_3(i,j))^2 / (2*Rreg*(TH_N_3(1,j)-TC_N_3(i,j)));
        else
            Qloss_reg_N_3(i,j)=(Troom_N_3-TC_N_3(i,j))/Rreg;
        end
    Qloss_tot_N_3(i,j)=Qloss_tubes_N_3(i,j)+Qloss_flowhead_N_3(i,j)+Qloss_reg_N_3(i,j)+Qloss_valve_N_3(1,j)+Qloss_filter_N_3(i,j);  %total heat losses
    QCloss_tot_N_3(i,j) = QC_N_3(i,j) - Qloss_tot_N_3(i,j);  %Cooling capacity simulated minus the losses
    end
end


%% Calculate losses from Numerical results for volume flow rate variation
%Old simulations:
Qloss_tubes_Cu_V=zeros(length(TC_V),length(V_V)); %K
Qloss_tubes_plas_V=zeros(length(TC_V),length(V_V)); %K
Qloss_flowhead_V=zeros(length(TC_V),length(V_V)); %K
Qloss_filter_V=zeros(length(TC_V),length(V_V)); %K
Qloss_valve_V=zeros(1,length(V_V)); %K
Qloss_reg_V=zeros(length(TC_V),length(V_V)); %K
Qloss_tot_V=zeros(length(TC_V),length(V_V)); %K

for i=1:length(TC_V)
    for j=1:length(V_V)
        Qloss_tubes_Cu_V(i,j)=(Troom_N-TC_V(i))/Rtubes_Cu_av; %W
        Qloss_tubes_plas_V(i,j)=(Troom_N-TC_V(i))/Rtubes_plas_av; %W
        Qloss_tubes_V(i,j)=Qloss_tubes_Cu_V(i,j) + Qloss_tubes_plas_V(i,j); %W
        Qloss_flowhead_V(i,j)=(Troom_N-TC_V(i))/Rflowhead; %W
        Qloss_filter_V(i,j)=(Troom_N-TC_V(i))/Rfilter_av; %W
        
%        Qloss_valve_V(1,j)=Tvalve_av*freq_V(:,j);  %CHECK! DO IT AGAIN!
        Qloss_valve_V(1,j)=13.2*freq_V(:,j)-2.5;
%         [valV,indTCV] = min(abs(TCout(:)-TC_V(i)));
%         Qvalves_V(1,j)=Tvalves(indTCV).*freq_V(:,j);
         
 
        if TC_V(i)<Troom_N
            Qloss_reg_V(i,j)=(Troom_N-TC_V(1,i))^2 / (2*Rreg*(TH_V-TC_V(i)));
        else
            Qloss_reg_V(i,j)=(Troom_N-TC_V(1,i))/Rreg;
        end        
      
        Qloss_tot_V(i,j)=Qloss_tubes_V(i,j)+Qloss_flowhead_V(i,j)+Qloss_reg_V(i,j)+Qloss_valve_V(1,j)+Qloss_filter_V(i,j);
        QCloss_tot_V(i,j) = QC_V(i,j) - Qloss_tot_V(i,j);
    end
end

%After changing the flow rate profile:
Qloss_tubes_Cu_V_2=zeros(length(TC_V_2),length(V_V_2)); %K
Qloss_tubes_plas_V_2=zeros(length(TC_V_2),length(V_V_2)); %K
Qloss_flowhead_V_2=zeros(length(TC_V_2),length(V_V_2)); %K
Qloss_filter_V_2=zeros(length(TC_V_2),length(V_V_2)); %K
Qloss_valve_V_2=zeros(1,length(V_V_2)); %K
Qloss_reg_V_2=zeros(length(TC_V_2),length(V_V_2)); %K
Qloss_tot_V_2=zeros(length(TC_V_2),length(V_V_2)); %K

for i=1:length(TC_V_2)
    for j=1:length(V_V_2)
        Qloss_tubes_Cu_V_2(i,j)=(Troom_N-TC_V_2(i))/Rtubes_Cu_av; %W
        Qloss_tubes_plas_V_2(i,j)=(Troom_N-TC_V_2(i))/Rtubes_plas_av; %W
        Qloss_tubes_V_2(i,j)=Qloss_tubes_Cu_V_2(i,j) + Qloss_tubes_plas_V_2(i,j); %W
        Qloss_flowhead_V_2(i,j)=(Troom_N-TC_V_2(i))/Rflowhead; %W
        Qloss_filter_V_2(i,j)=(Troom_N-TC_V_2(i))/Rfilter_av; %W
        
%        Qloss_valve_V_2(1,j)=Tvalve_av*freq_V_2(:,j);  %CHECK! DO IT AGAIN!
        Qloss_valve_V_2(1,j)=13.2*freq_V_2(:,j)-5;
%         [valV,indTCV] = min(abs(TCout(:)-TC_V_2(i)));
%         Qvalves_V_2(1,j)=Tvalves(indTCV).*freq_V_2(:,j);
     
        if TC_V_2(i)<Troom_N
            Qloss_reg_V_2(i,j)=(Troom_N-TC_V_2(1,i))^2 / (2*Rreg*(TH_V_2-TC_V_2(i)));
        else
            Qloss_reg_V_2(i,j)=(Troom_N-TC_V_2(1,i))/Rreg;
        end        
      
        Qloss_tot_V_2(i,j)=Qloss_tubes_V_2(i,j)+Qloss_flowhead_V_2(i,j)+Qloss_reg_V_2(i,j)+Qloss_valve_V_2(1,j)+Qloss_filter_V_2(i,j);
        QCloss_tot_V_2(i,j) = QC_V_2(i,j) - Qloss_tot_V_2(i,j);
    end
end

%% Calculate losses from Numerical results for volume flow rate variation and different particle size:

%After changing the flow rate profile for 0.45 mm particle size:
Qloss_tubes_Cu_V_045=zeros(length(TC_V_045),length(V_V_045)); %K
Qloss_tubes_plas_V_045=zeros(length(TC_V_045),length(V_V_045)); %K
Qloss_flowhead_V_045=zeros(length(TC_V_045),length(V_V_045)); %K
Qloss_filter_V_045=zeros(length(TC_V_045),length(V_V_045)); %K
Qloss_valve_V_045=zeros(1,length(V_V_045)); %K
Qloss_reg_V_045=zeros(length(TC_V_045),length(V_V_045)); %K
Qloss_tot_V_045=zeros(length(TC_V_045),length(V_V_045)); %K

for i=1:length(TC_V_045)
    for j=1:length(V_V_045)
        Qloss_tubes_Cu_V_045(i,j)=(Troom_N-TC_V_045(i))/Rtubes_Cu_av; %W
        Qloss_tubes_plas_V_045(i,j)=(Troom_N-TC_V_045(i))/Rtubes_plas_av; %W
        Qloss_tubes_V_045(i,j)=Qloss_tubes_Cu_V_045(i,j) + Qloss_tubes_plas_V_045(i,j); %W
        Qloss_flowhead_V_045(i,j)=(Troom_N-TC_V_045(i))/Rflowhead; %W
        Qloss_filter_V_045(i,j)=(Troom_N-TC_V_045(i))/Rfilter_av; %W
        
%        Qloss_valve_V_045(1,j)=Tvalve_av*freq_V_045(:,j);  %CHECK! DO IT AGAIN!
        Qloss_valve_V_045(1,j)=13.2*freq_V_045(:,j)-5;
%         [valV,indTCV] = min(abs(TCout(:)-TC_V_045(i)));
%         Qvalves_V_045(1,j)=Tvalves(indTCV).*freq_V_045(:,j);
     
        if TC_V_045(i)<Troom_N
            Qloss_reg_V_045(i,j)=(Troom_N-TC_V_045(1,i))^2 / (2*Rreg*(TH_V_045-TC_V_045(i)));
        else
            Qloss_reg_V_045(i,j)=(Troom_N-TC_V_045(1,i))/Rreg;
        end        
      
        Qloss_tot_V_045(i,j)=Qloss_tubes_V_045(i,j)+Qloss_flowhead_V_045(i,j)+Qloss_reg_V_045(i,j)+Qloss_valve_V_045(1,j)+Qloss_filter_V_045(i,j);
        QCloss_tot_V_045(i,j) = QC_V_045(i,j) - Qloss_tot_V_045(i,j);
    end
end

dTloss_int_V_045 = zeros(6,length(V_V_045));
for i=1:length(V_V_045)
    dTloss_int_V_045(1,i) = interp1(QCloss_tot_V_045(:,i),dT_V_045(:),0);
    dTloss_int_V_045(2,i) = interp1(QCloss_tot_V_045(:,i),dT_V_045(:),100);
    dTloss_int_V_045(3,i) = interp1(QCloss_tot_V_045(:,i),dT_V_045(:),200);
    dTloss_int_V_045(4,i) = interp1(QCloss_tot_V_045(:,i),dT_V_045(:),300);
    dTloss_int_V_045(5,i) = interp1(QCloss_tot_V_045(:,i),dT_V_045(:),400);
    dTloss_int_V_045(6,i) = interp1(QCloss_tot_V_045(:,i),dT_V_045(:),500);
end

%After changing the flow rate profile for 0.5 mm particle size:
Qloss_tubes_Cu_V_05=zeros(length(TC_V_05),length(V_V_05)); %K
Qloss_tubes_plas_V_05=zeros(length(TC_V_05),length(V_V_05)); %K
Qloss_flowhead_V_05=zeros(length(TC_V_05),length(V_V_05)); %K
Qloss_filter_V_05=zeros(length(TC_V_05),length(V_V_05)); %K
Qloss_valve_V_05=zeros(1,length(V_V_05)); %K
Qloss_reg_V_05=zeros(length(TC_V_05),length(V_V_05)); %K
Qloss_tot_V_05=zeros(length(TC_V_05),length(V_V_05)); %K

for i=1:length(TC_V_05)
    for j=1:length(V_V_05)
        Qloss_tubes_Cu_V_05(i,j)=(Troom_N-TC_V_05(i))/Rtubes_Cu_av; %W
        Qloss_tubes_plas_V_05(i,j)=(Troom_N-TC_V_05(i))/Rtubes_plas_av; %W
        Qloss_tubes_V_05(i,j)=Qloss_tubes_Cu_V_05(i,j) + Qloss_tubes_plas_V_05(i,j); %W
        Qloss_flowhead_V_05(i,j)=(Troom_N-TC_V_05(i))/Rflowhead; %W
        Qloss_filter_V_05(i,j)=(Troom_N-TC_V_05(i))/Rfilter_av; %W
        
%        Qloss_valve_V_05(1,j)=Tvalve_av*freq_V_05(:,j);  %CHECK! DO IT AGAIN!
        Qloss_valve_V_05(1,j)=13.2*freq_V_05(:,j)-5;
%         [valV,indTCV] = min(abs(TCout(:)-TC_V_05(i)));
%         Qvalves_V_05(1,j)=Tvalves(indTCV).*freq_V_05(:,j);
     
        if TC_V_05(i)<Troom_N
            Qloss_reg_V_05(i,j)=(Troom_N-TC_V_05(1,i))^2 / (2*Rreg*(TH_V_05-TC_V_05(i)));
        else
            Qloss_reg_V_05(i,j)=(Troom_N-TC_V_05(1,i))/Rreg;
        end        
      
        Qloss_tot_V_05(i,j)=Qloss_tubes_V_05(i,j)+Qloss_flowhead_V_05(i,j)+Qloss_reg_V_05(i,j)+Qloss_valve_V_05(1,j)+Qloss_filter_V_05(i,j);
        QCloss_tot_V_05(i,j) = QC_V_05(i,j) - Qloss_tot_V_05(i,j);
    end
end

dTloss_int_V_05 = zeros(6,length(V_V_05));
for i=1:length(V_V_05)
    dTloss_int_V_05(1,i) = interp1(QCloss_tot_V_05(:,i),dT_V_05(:),0);
    dTloss_int_V_05(2,i) = interp1(QCloss_tot_V_05(:,i),dT_V_05(:),100);
    dTloss_int_V_05(3,i) = interp1(QCloss_tot_V_05(:,i),dT_V_05(:),200);
    dTloss_int_V_05(4,i) = interp1(QCloss_tot_V_05(:,i),dT_V_05(:),300);
    dTloss_int_V_05(5,i) = interp1(QCloss_tot_V_05(:,i),dT_V_05(:),400);
    dTloss_int_V_05(6,i) = interp1(QCloss_tot_V_05(:,i),dT_V_05(:),500);
end

%After changing the flow rate profile for 0.7 mm particle size:
Qloss_tubes_Cu_V_07=zeros(length(TC_V_07),length(V_V_07)); %K
Qloss_tubes_plas_V_07=zeros(length(TC_V_07),length(V_V_07)); %K
Qloss_flowhead_V_07=zeros(length(TC_V_07),length(V_V_07)); %K
Qloss_filter_V_07=zeros(length(TC_V_07),length(V_V_07)); %K
Qloss_valve_V_07=zeros(1,length(V_V_07)); %K
Qloss_reg_V_07=zeros(length(TC_V_07),length(V_V_07)); %K
Qloss_tot_V_07=zeros(length(TC_V_07),length(V_V_07)); %K

for i=1:length(TC_V_07)
    for j=1:length(V_V_07)
        Qloss_tubes_Cu_V_07(i,j)=(Troom_N-TC_V_07(i))/Rtubes_Cu_av; %W
        Qloss_tubes_plas_V_07(i,j)=(Troom_N-TC_V_07(i))/Rtubes_plas_av; %W
        Qloss_tubes_V_07(i,j)=Qloss_tubes_Cu_V_07(i,j) + Qloss_tubes_plas_V_07(i,j); %W
        Qloss_flowhead_V_07(i,j)=(Troom_N-TC_V_07(i))/Rflowhead; %W
        Qloss_filter_V_07(i,j)=(Troom_N-TC_V_07(i))/Rfilter_av; %W
        
%        Qloss_valve_V_07(1,j)=Tvalve_av*freq_V_07(:,j);  %CHECK! DO IT AGAIN!
        Qloss_valve_V_07(1,j)=13.2*freq_V_07(:,j)-5;
%         [valV,indTCV] = min(abs(TCout(:)-TC_V_07(i)));
%         Qvalves_V_07(1,j)=Tvalves(indTCV).*freq_V_07(:,j);
     
        if TC_V_07(i)<Troom_N
            Qloss_reg_V_07(i,j)=(Troom_N-TC_V_07(1,i))^2 / (2*Rreg*(TH_V_07-TC_V_07(i)));
        else
            Qloss_reg_V_07(i,j)=(Troom_N-TC_V_07(1,i))/Rreg;
        end        
      
        Qloss_tot_V_07(i,j)=Qloss_tubes_V_07(i,j)+Qloss_flowhead_V_07(i,j)+Qloss_reg_V_07(i,j)+Qloss_valve_V_07(1,j)+Qloss_filter_V_07(i,j);
        QCloss_tot_V_07(i,j) = QC_V_07(i,j) - Qloss_tot_V_07(i,j);
    end
end

dTloss_int_V_07 = zeros(6,length(V_V_07));
for i=1:length(V_V_07)
    dTloss_int_V_07(1,i) = interp1(QCloss_tot_V_07(:,i),dT_V_07(:),0);
    dTloss_int_V_07(2,i) = interp1(QCloss_tot_V_07(:,i),dT_V_07(:),100);
    dTloss_int_V_07(3,i) = interp1(QCloss_tot_V_07(:,i),dT_V_07(:),200);
    dTloss_int_V_07(4,i) = interp1(QCloss_tot_V_07(:,i),dT_V_07(:),300);
    dTloss_int_V_07(5,i) = interp1(QCloss_tot_V_07(:,i),dT_V_07(:),400);
    dTloss_int_V_07(6,i) = interp1(QCloss_tot_V_07(:,i),dT_V_07(:),500);
end

%After changing the flow rate profile for 0.75 mm particle size:
Qloss_tubes_Cu_V_075=zeros(length(TC_V_075),length(V_V_075)); %K
Qloss_tubes_plas_V_075=zeros(length(TC_V_075),length(V_V_075)); %K
Qloss_flowhead_V_075=zeros(length(TC_V_075),length(V_V_075)); %K
Qloss_filter_V_075=zeros(length(TC_V_075),length(V_V_075)); %K
Qloss_valve_V_075=zeros(1,length(V_V_075)); %K
Qloss_reg_V_075=zeros(length(TC_V_075),length(V_V_075)); %K
Qloss_tot_V_075=zeros(length(TC_V_075),length(V_V_075)); %K

for i=1:length(TC_V_075)
    for j=1:length(V_V_075)
        Qloss_tubes_Cu_V_075(i,j)=(Troom_N-TC_V_075(i))/Rtubes_Cu_av; %W
        Qloss_tubes_plas_V_075(i,j)=(Troom_N-TC_V_075(i))/Rtubes_plas_av; %W
        Qloss_tubes_V_075(i,j)=Qloss_tubes_Cu_V_075(i,j) + Qloss_tubes_plas_V_075(i,j); %W
        Qloss_flowhead_V_075(i,j)=(Troom_N-TC_V_075(i))/Rflowhead; %W
        Qloss_filter_V_075(i,j)=(Troom_N-TC_V_075(i))/Rfilter_av; %W
        
%        Qloss_valve_V_075(1,j)=Tvalve_av*freq_V_075(:,j);  %CHECK! DO IT AGAIN!
        Qloss_valve_V_075(1,j)=13.2*freq_V_075(:,j)-5;
%         [valV,indTCV] = min(abs(TCout(:)-TC_V_075(i)));
%         Qvalves_V_075(1,j)=Tvalves(indTCV).*freq_V_075(:,j);
     
        if TC_V_075(i)<Troom_N
            Qloss_reg_V_075(i,j)=(Troom_N-TC_V_075(1,i))^2 / (2*Rreg*(TH_V_075-TC_V_075(i)));
        else
            Qloss_reg_V_075(i,j)=(Troom_N-TC_V_075(1,i))/Rreg;
        end        
      
        Qloss_tot_V_075(i,j)=Qloss_tubes_V_075(i,j)+Qloss_flowhead_V_075(i,j)+Qloss_reg_V_075(i,j)+Qloss_valve_V_075(1,j)+Qloss_filter_V_075(i,j);
        QCloss_tot_V_075(i,j) = QC_V_075(i,j) - Qloss_tot_V_075(i,j);
    end
end

dTloss_int_V_075 = zeros(6,length(V_V_075));
for i=1:length(V_V_075)
    dTloss_int_V_075(1,i) = interp1(QCloss_tot_V_075(:,i),dT_V_075(:),0);
    dTloss_int_V_075(2,i) = interp1(QCloss_tot_V_075(:,i),dT_V_075(:),100);
    dTloss_int_V_075(3,i) = interp1(QCloss_tot_V_075(:,i),dT_V_075(:),200);
    dTloss_int_V_075(4,i) = interp1(QCloss_tot_V_075(:,i),dT_V_075(:),300);
    dTloss_int_V_075(5,i) = interp1(QCloss_tot_V_075(:,i),dT_V_075(:),400);
    dTloss_int_V_075(6,i) = interp1(QCloss_tot_V_075(:,i),dT_V_075(:),500);
end


%% Calculate losses from Numerical results for frequency variation

Qloss_tubes_Cu_F=zeros(length(TC_F),length(freq_F)); %K
Qloss_tubes_plas_F=zeros(length(TC_F),length(freq_F)); %K
Qloss_flowhead_F=zeros(length(TC_F),length(freq_F)); %K
Qloss_filter_F=zeros(length(TC_F),length(freq_F)); %K
Qloss_valve_F=zeros(1,length(freq_F)); %K
Qloss_reg_F=zeros(length(TC_F),length(freq_F)); %K
Qloss_tot_F=zeros(length(TC_F),length(freq_F)); %K

for i=1:length(TC_F)
    for j=1:length(freq_F)
        Qloss_tubes_Cu_F(i,j)=(Troom_N-TC_F(i))/Rtubes_Cu_av; %W
        Qloss_tubes_plas_F(i,j)=(Troom_N-TC_F(i))/Rtubes_plas_av; %W
        Qloss_tubes_F(i,j)=Qloss_tubes_Cu_F(i,j) + Qloss_tubes_plas_F(i,j); %W
        Qloss_flowhead_F(i,j)=(Troom_N-TC_F(i))/Rflowhead; %W
        Qloss_filter_F(i,j)=(Troom_N-TC_F(i))/Rfilter_av; %W

        
        Qloss_valve_F(1,j)=Tvalve*freq_F(:,j);  %CHECK! DO IT AGAIN!
        %Qloss_valve_F(1,j)=13.2*freq_F(:,j)-2.5;
%         [valV,indTCV] = min(abs(TCout(:)-TC_V(i)));
%         Qvalves_V(1,j)=Tvalves(indTCV).*freq_V(:,j);
         
      
        if TC_F(i)<Troom_N
            Qloss_reg_F(i,j)=(Troom_N-TC_F(1,i))^2 / (2*Rreg*(TH_F-TC_F(i)));
        else
            Qloss_reg_F(i,j)=(Troom_N-TC_F(1,i))/Rreg;
        end        
      
        Qloss_tot_F(i,j)=Qloss_tubes_F(i,j)+Qloss_flowhead_F(i,j)+Qloss_reg_F(i,j)+Qloss_valve_F(1,j)+Qloss_filter_F(i,j);
        QCloss_tot_F(i,j) = QC_F(i,j) - Qloss_tot_F(i,j);
    end
end

%% Find Temperature span at certain Cooling Capacity from Numerical results

dTloss_int_N = zeros(6,length(TH_N));
QCloss_tot_N(1,2)=600;  % This interpolation has the problem that the data is not descendent (continue). So it needs to be fix by interpolating only with decreasing cooling capacity. 
QCloss_tot_N(1,1)=480;
QCloss_tot_N(2,1)=450;
QCloss_tot_N(1,3)=600;
QCloss_tot_N(2,3)=560;
for i=1:length(TH_N)
    dTloss_int_N(1,i) = interp1(QCloss_tot_N(:,i),dT_N(:),0);
    dTloss_int_N(2,i) = interp1(QCloss_tot_N(:,i),dT_N(:),100);
    dTloss_int_N(3,i) = interp1(QCloss_tot_N(:,i),dT_N(:),200);
    dTloss_int_N(4,i) = interp1(QCloss_tot_N(:,i),dT_N(:),300);
    dTloss_int_N(5,i) = interp1(QCloss_tot_N(:,i),dT_N(:),400);
    dTloss_int_N(6,i) = interp1(QCloss_tot_N(:,i),dT_N(:),500);
end
dTloss_int_N(2,10)=nan;

%After changing mass flow rate profile (400L/h):
dTloss_int_N_2 = zeros(6,length(TH_N_2));
%QCloss_tot_N_2(1,2)=600;  % This interpolation has the problem that the data is not descendent (continue). So it needs to be fix by interpolating only with decreasing cooling capacity. 
%QCloss_tot_N_2(1,1)=480;
%QCloss_tot_N_2(2,1)=450;
%QCloss_tot_N_2(1,3)=600;
%QCloss_tot_N_2(2,3)=560;
for i=1:length(TH_N_2)
    dTloss_int_N_2(1,i) = interp1(QCloss_tot_N_2(:,i),dT_N_2(:),0);
    dTloss_int_N_2(2,i) = interp1(QCloss_tot_N_2(:,i),dT_N_2(:),100);
    dTloss_int_N_2(3,i) = interp1(QCloss_tot_N_2(:,i),dT_N_2(:),200);
    dTloss_int_N_2(4,i) = interp1(QCloss_tot_N_2(:,i),dT_N_2(:),300);
    %dTloss_int_N_2(5,i) = interp1(QCloss_tot_N_2(:,i),dT_N_2(:),400,'linear','extrap');
    dTloss_int_N_2(5,i) = interp1(QCloss_tot_N_2(:,i),dT_N_2(:),400);
    dTloss_int_N_2(6,i) = interp1(QCloss_tot_N_2(:,i),dT_N_2(:),500);
end
%dTloss_int_N_2(2,10)=nan;
%dTloss_int_N_2(5,1) = 1.5;

%After changing mass flow rate profile (440L/h):
dTloss_int_N_3 = zeros(6,length(TH_N_3));
%QCloss_tot_N_3(1,2)=600;  % This interpolation has the problem that the data is not descendent (continue). So it needs to be fix by interpolating only with decreasing cooling capacity. 
%QCloss_tot_N_3(1,1)=480;
%QCloss_tot_N_3(2,1)=450;
%QCloss_tot_N_3(1,3)=600;
%QCloss_tot_N_3(2,3)=560;
for i=1:length(TH_N_3)
    dTloss_int_N_3(1,i) = interp1(QCloss_tot_N_3(:,i),dT_N_3(:),0);
    dTloss_int_N_3(2,i) = interp1(QCloss_tot_N_3(:,i),dT_N_3(:),100);
    dTloss_int_N_3(3,i) = interp1(QCloss_tot_N_3(:,i),dT_N_3(:),200);
    dTloss_int_N_3(4,i) = interp1(QCloss_tot_N_3(:,i),dT_N_3(:),300);
    dTloss_int_N_3(5,i) = interp1(QCloss_tot_N_3(:,i),dT_N_3(:),400);
    dTloss_int_N_3(6,i) = interp1(QCloss_tot_N_3(:,i),dT_N_3(:),500);
end
%dTloss_int_N_3(2,10)=nan;

dTloss_int_V = zeros(6,length(V_V));
for i=1:length(V_V)
    dTloss_int_V(1,i) = interp1(QCloss_tot_V(:,i),dT_V(:),0);
    dTloss_int_V(2,i) = interp1(QCloss_tot_V(:,i),dT_V(:),100);
    dTloss_int_V(3,i) = interp1(QCloss_tot_V(:,i),dT_V(:),200);
    dTloss_int_V(4,i) = interp1(QCloss_tot_V(:,i),dT_V(:),300);
    dTloss_int_V(5,i) = interp1(QCloss_tot_V(:,i),dT_V(:),400);
    dTloss_int_V(6,i) = interp1(QCloss_tot_V(:,i),dT_V(:),500);
end
%After changing the flow profile:
dTloss_int_V_2 = zeros(6,length(V_V_2));
for i=1:length(V_V_2)
    dTloss_int_V_2(1,i) = interp1(QCloss_tot_V_2(:,i),dT_V_2(:),0);
    dTloss_int_V_2(2,i) = interp1(QCloss_tot_V_2(:,i),dT_V_2(:),100);
    dTloss_int_V_2(3,i) = interp1(QCloss_tot_V_2(:,i),dT_V_2(:),200);
    dTloss_int_V_2(4,i) = interp1(QCloss_tot_V_2(:,i),dT_V_2(:),300);
    dTloss_int_V_2(5,i) = interp1(QCloss_tot_V_2(:,i),dT_V_2(:),400);
    dTloss_int_V_2(6,i) = interp1(QCloss_tot_V_2(:,i),dT_V_2(:),500);
end
dTloss_int_V_2(5,6) = -2.42; %extrapolation from QCloss results - aprox. 1 K difference due to about 10 W loss (flow distributor + gain cold end)
%dTloss_int_V_2(5,6) = -1.05; %extrapolation from dTloss results
%dTloss_int_V_2(5,6) = (-2.42 - 1.05) / 2;

dTloss_int_F = zeros(6,length(freq_F));
for i=1:length(freq_F)
    dTloss_int_F(1,i) = interp1(QCloss_tot_F(:,i),dT_F(:),0);
    dTloss_int_F(2,i) = interp1(QCloss_tot_F(:,i),dT_F(:),100);
    dTloss_int_F(3,i) = interp1(QCloss_tot_F(:,i),dT_F(:),200);
    dTloss_int_F(4,i) = interp1(QCloss_tot_F(:,i),dT_F(:),300);
    dTloss_int_F(5,i) = interp1(QCloss_tot_F(:,i),dT_F(:),400);
    dTloss_int_F(6,i) = interp1(QCloss_tot_F(:,i),dT_F(:),500);
end

%% Calculate loss in VALVE:
%For hot in temperature dependence:
dTloss_int_N_novalve = zeros(2,length(TH_N));
%deltaT_valve = zeros(5,13);
for i=1:length(TC_N)
    for j=1:length(TH_N)
        Qloss_tot_N_novalve(i,j)=Qloss_tubes_N(i,j)+Qloss_flowhead_N(i,j)+Qloss_reg_N(i,j)+Qloss_filter_N(i,j);  %total heat losses
        QCloss_tot_N_novalve(i,j) = QC_N(i,j) - Qloss_tot_N_novalve(i,j);  %Cooling capacity simulated minus the losses (without loss in valve)
    end
end

QCloss_tot_N_novalve(1,1)=490;
QCloss_tot_N_novalve(2,1)=460;
QCloss_tot_N_novalve(1,2)=600;
QCloss_tot_N_novalve(2,2)=560;
QCloss_tot_N_novalve(1,3)=620;
QCloss_tot_N_novalve(2,3)=570;

for i=1:length(TH_N)
    dTloss_int_N_novalve(1,i) = interp1(QCloss_tot_N_novalve(:,i),dT_N(:),200);
    dTloss_int_N_novalve(2,i) = interp1(QCloss_tot_N_novalve(:,i),dT_N(:),400);
end
    
deltaT = THout-TCout;

for i=1:11
    if i < 6
        deltaT_valve(1,i) = dTloss_int_N_novalve(1,i) - deltaT(1,12-i);
    else 
        deltaT_valve(1,i) = dTloss_int_N_novalve(1,i+1) - deltaT(1,12-i);
    end
end

for i=1:9
    deltaT_valve(2,i) = dTloss_int_N_novalve(2,i) - deltaT(2,i);
end

for i=1:11
    Qloss_valve_N_dT(1,i) = deltaT_valve(1,i) * 200 / deltaT(1,12-i);
    Qloss_valve_N_dT(2,i) = deltaT_valve(2,i) * 400 / deltaT(2,i);
end    

%Calculation of Tvalve dependence on TC:
%In Excel: TCout2(1,:),Qloss_valve_N_dT(1,:)

% For volume flow rate dependence:
dTloss_int_V_novalve = zeros(2,length(TH_N));
for i=1:length(TC_V)
    for j=1:length(V_V)
        Qloss_tot_V_novalve(i,j)=Qloss_tubes_V(i,j)+Qloss_flowhead_V(i,j)+Qloss_reg_V(i,j)+Qloss_filter_V(i,j);
        QCloss_tot_V_novalve(i,j) = QC_V(i,j) - Qloss_tot_V_novalve(i,j);        
    end
end
dTloss_int_V_novalve = zeros(2,length(V_V));
for i=1:length(V_V)
    dTloss_int_V_novalve(1,i) = interp1(QCloss_tot_V_novalve(:,i),dT_V(:),200);
    dTloss_int_V_novalve(2,i) = interp1(QCloss_tot_V_novalve(:,i),dT_V(:),400);
end

for i=1:5
    if i<4
        deltaT_valve(3,i) = dTloss_int_V_novalve(1,i+1) - deltaT(3,i);
    else
        deltaT_valve(3,i) = dTloss_int_V_novalve(1,i+2) - deltaT(3,i);
    end
end

for i=1:4
    deltaT_valve(4,i) = dTloss_int_V_novalve(2,i+4) - deltaT(4,5-i);
end

deltaT_valve(deltaT_valve==0) = nan; %plot nonzero values

for i=1:5
    Qloss_valve_V_dT(1,i) = deltaT_valve(3,i) * 200 / deltaT(3,i);
end    
for i=1:4
    Qloss_valve_V_dT(2,i) = deltaT_valve(4,i) * 400 / deltaT(4,5-i);
end    

%For frequency dependence:
dTloss_int_F_novalve = zeros(1,length(freq_F));
for i=1:length(TC_F)
    for j=1:length(freq_F)
        Qloss_tot_F_novalve(i,j)=Qloss_tubes_F(i,j)+Qloss_flowhead_F(i,j)+Qloss_reg_F(i,j)+Qloss_filter_F(i,j);  %total heat losses
        QCloss_tot_F_novalve(i,j) = QC_F(i,j) - Qloss_tot_F_novalve(i,j);  %Cooling capacity simulated minus the losses (without loss in valve)
    end
end

for i=1:length(freq_F)
    dTloss_int_F_novalve(1,i) = interp1(QCloss_tot_F_novalve(:,i),dT_F(:),200);
    TC_F_novalve(1,i) = TH_F - dTloss_int_F_novalve(1,i);
end

deltaT_valve(5,1) = dTloss_int_F_novalve(1,1) - deltaT(5,2);
Qloss_valve_F_dT(1,1) = deltaT_valve(5,1) * 200 / deltaT(5,2);
deltaT_valve(5,2) = dTloss_int_F_novalve(1,2) - deltaT(5,3);
Qloss_valve_F_dT(1,2) = deltaT_valve(5,2) * 200 / deltaT(5,3);
deltaT_valve(5,3) = dTloss_int_F_novalve(1,3) - deltaT(5,5);
Qloss_valve_F_dT(1,3) = deltaT_valve(5,3) * 200 / deltaT(5,5);
deltaT_valve(5,4) = dTloss_int_F_novalve(1,4) - deltaT(5,7);
Qloss_valve_F_dT(1,4) = deltaT_valve(5,4) * 200 / deltaT(5,7);
deltaT_valve(5,5) = dTloss_int_F_novalve(1,5) - deltaT(5,9);
Qloss_valve_F_dT(1,5) = deltaT_valve(5,5) * 200 / deltaT(5,9);
deltaT_valve(5,6) = dTloss_int_F_novalve(1,6) - deltaT(5,11);
Qloss_valve_F_dT(1,6) = deltaT_valve(5,6) * 200 / deltaT(5,11);
deltaT_valve(5,7) = dTloss_int_F_novalve(1,7) - deltaT(5,13);
Qloss_valve_F_dT(1,7) = deltaT_valve(5,7) * 200 / deltaT(5,13);

%plot(THin(1,1:11),Qloss_reg(1,1:11))
%Calculation of Tvalve dependence on TC and f:

%In Excel: TCout2(1,:),Qloss_valve_N_dT(1,:)

%% Performance Metrics (experimental data)

DeltaT = THout-TCout;
Deltap = pHin-pHout;
eta_pump = 0.7;
Wpump = V.* Deltap./(36*eta_pump); %V in L/h and DeltaP in bar
W = Wmotor + Wpump;  %Wmotor is experimental, Wpump is calculated
for i=1:length(Wpump(:,1))
    Wpump(Wpump==0) = NaN;
    Wpump_av(i,:) = nanmean(Wpump(i,:)); %Average pumping power for set of experiments
end

COP = QC./W; %Experimental COP
COP_noloss = (QC + Qloss_tot_exp) ./ W; %COP experimental without thermal losses 
%COP_novalve = (QC + Qloss_valve*2) ./ W; %COP experimental without friction at the valves

%COP analysis taking into account the efficiency of the motor: 
W_system = 220; % W - Power of the motor at 1.5 Hz for the whole system
W_noload = 133; % W - Power of the motor at 1.5 Hz for the motor spining only with the frequency inversor
W_novalves = 154; % W - Power of the motor at 1.5 Hz for the motor spining only with the AMR but without the valves on

%W_AMR = W_system - W_noload; %at 1.5 Hz - Work that is needed to rotate the whole AMR (magnetic work) + the flowhead inside the valves (friction)
Wmotor(Wmotor==0)= nan;
W_AMR(1:4,:) = Wmotor(1:4,:) - W_noload; %at 1.5 Hz - Work that is needed to rotate the whole AMR (magnetic work) + the flowhead inside the valves (friction) - from experimental results of the power of the motor
Eff_motor = 0.8; %Electric motor nominal efficiency
% The other COP to plot:
COP_novalve(1:4,:) = QC(1:4,:) ./ (Wpump(1:4,:) + Wmotor(1:4,:) - (W_system - W_novalves)) ; %COP experimental without friction at the valves at 1.5 Hz
COP_AMR(1:4,:) = QC(1:4,:) ./ (Wpump(1:4,:) + W_AMR(1:4,:));
COP_motor_eff(1:4,:) = QC(1:4,:) ./ (Wpump(1:4,:) + (W_AMR(1:4,:)/Eff_motor)) ; %COP experimental without friction at the valves at 1.5 Hz

W_mag = W_novalves - W_noload; % W - Power of the motor at 1.5 Hz for the magnetic work (including eddy currents) 
%COP_cycle(1:4,:) = QC(1:4,:) ./ (Wpump(1:4,:) + W_mag);
COP_cycle(1:4,:) = QC(1:4,:) ./ ((Wpump(1:4,:)*0.7) + (W_mag*0.8));

Gain_COP_novalve(1:4,:) = (COP_novalve(1:4,:)- COP(1:4,:)) ./ COP_novalve(1:4,:);
Gain_COP_novalve_av = nanmean(nanmean(Gain_COP_novalve)); 
Gain_COP_motor_eff(1:4,:) = (COP_motor_eff(1:4,:)- COP(1:4,:)) ./ COP_motor_eff(1:4,:);
Gain_COP_motor_av = nanmean(nanmean(Gain_COP_motor_eff)); 
Gain_COP_AMR(1:4,:) = (COP_AMR(1:4,:)- COP(1:4,:)) ./ COP_AMR(1:4,:);
Gain_COP_AMR_av = nanmean(nanmean(Gain_COP_AMR)); 

for i=1:length(COP(:,1))
    COP_av(i,:) = nanmean(COP(i,:)); %Average pumping power for set of experiments
end

%Calculate COP Carnot by the temperatures of the fluid:
rho_f = 1028; %kg/m3
Cp_f = 3903; %kg/JK
QC_calc = V.*(TCin-TCout)*rho_f*Cp_f/(3600*1000); %W
QH_calc = V.*(THout-THin)*rho_f*Cp_f/(3600*1000); %W
COP_Qcalc = QC_calc./(QH_calc-QC_calc);
COP_QCcalc = QC./(QH_calc-QC);

% Calculate LMTD for 'ambient' temperature: 
DeltaTH_LMTD = ((THout-TChiller)-(THin-TChiller))./log((THout-TChiller)./(THin-TChiller));

COP_ideal_e = TCout ./ ((TChiller + DeltaTH_LMTD) - TCout);
COP_ideal = TCout ./ (THout - TCout);

% COP_carnot_Tch = TCout ./ (TChiller  - TCout)
% COP_carnot_TH = TCout ./ (THout - TCout)

%Ex_Q = QC.*((Troom./TCout) - 1); % - Exergetic... - Rowe Performance Metrics AMR
Ex_Q = QC.*((THout./TCout) - 1); % - Exergetic... - Russek and Arnold
eta_Ex = Ex_Q./W; %Efficiency defined in terms of Exergy - Rowe Performance Metrics AMR 
%mu_Ex200(1,:) = Ex_Q200(1,:)./(B_0*V_MCM) %Specific Exergetic Cooling Power - Arnold

% Split 2nd law efficiency in Internal and External (Jader)
% It can be separate the Heat Ex. losses (external irrev)
% this is the baseline to compare with other technologies

eta_system = COP./COP_ideal; %Carnot uses THout
eta_system_e = COP./COP_ideal_e; %Carnot uses THout
%eta_2ndlaw_Tch = COP./COP_carnot_Tch; %Carnot uses THout
eta_cycle = COP_cycle ./ COP_ideal(1:4,:); %Carnot uses THout

for i=1:length(eta_system(:,1))
    eta_system_av(i,:) = nanmean(eta_system(i,:))
end

%eta_i =  Paper Hermes & Jader 
%eta_e =  Paper Hermes & Jader

%what is the role of the internal and external irreversibilites
%for the COP and overall performance

%break up COP
%Exergy?

%% Calculate Performance metrics for numerical data:

COP_N = zeros(2,length(dTint));
for i=1:length(dTint(1,:))
    COP_N(1,i) = (TH_N(1,i) - dTint(1,i)) ./ (dTint(1,i));        
    COP_N(2,i) = (TH_N(1,i) - dTint(2,i)) ./ (dTint(2,i));
end

COP_V = zeros(2,length(dTint_V));
for i=1:length(dTint_V(1,:))
    COP_V(1,:) = (TH_V - dTint_V(1,:)) ./ (dTint_V(1,:));        
    COP_V(2,:) = (TH_V - dTint_V(2,:)) ./ (dTint_V(2,:));
end

%COP_V = TC_V ./ (TH_V-TC_V);

%% PLOTS Temperature span 
% %Hot end in temperature dependence
% % 
THin(THin==273) = nan; %plot nonzero values
%deltaT(deltaT==0) = nan; %plot nonzero values

% figure1 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure1,'FontSize',14,'FontName','Times New Roman');
% xlim([285 305]);
% ylim([0 25]);
% box('on');
% grid('on');
% hold('all');
% plot(THin(1,:),deltaT(1,:),'bo','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',7)
% plot(TH_N_2(:),dTint_N_2(1,:),':b','Linewidth',2','MarkerSize',6)
% plot(TH_N_2(:),dTloss_int_N_2(3,:),'-b','Linewidth',1.5','MarkerFaceColor','w','MarkerSize',5)
% % plot(TH_N(:),dTint(1,:),':bo','Linewidth',1.5','MarkerSize',6)
% % plot(TH_N(:),dTloss_int_N(3,:),'--bo','Linewidth',1.5','MarkerSize',6)
% %plot(TH_N(:),dTloss_int_N_novalve(1,:),'--co','Linewidth',1.5','MarkerSize',6)
% plot(THin(2,:),deltaT(2,:),'r^','Linewidth',2','MarkerFaceColor','r','MarkerSize',7)
% plot(TH_N_2(:),dTint_N_2(2,:),':r','Linewidth',1.5','MarkerSize',7)
% plot(TH_N_2(:),dTloss_int_N_2(5,:),'-r','Linewidth',1.5','MarkerFaceColor','w','MarkerSize',6)
% % plot(TH_N(:),dTint(2,:),':r^','Linewidth',1.5','MarkerSize',6)
% % plot(TH_N(1,1:9),dTloss_int_N(5,1:9),'--r^','Linewidth',1.5','MarkerSize',6)
% %plot(TH_N(1,1:9),dTloss_int_N_novalve(2,1:9),'--c^','Linewidth',1.5','MarkerSize',6)
% print('-depsc','fig_DeltaT.eps');
% xlabel('Hot reservoir temperature (K)', 'Fontsize', 16,'FontName','Times New Roman')
% ylabel('Temperature Span (K)', 'Fontsize', 16,'FontName','Times New Roman')
% %h1 = legend('$\dot{\mathrm{Q}}_{C} =$ 200 W (Experimental)','$\dot{\mathrm{Q}}_{C} =$ 200 W (Numerical)','$\dot{\mathrm{Q}}_{C} =$ 200 W (Num. Losses)','$\dot{\mathrm{Q}}_{C} =$ 400 W (Experimental)','$\dot{\mathrm{Q}}_{C} =$ 400 W (Numerical)','$\dot{\mathrm{Q}}_{C} =$ 400 W (Num. Losses)');
% h1 = legend('200 W (Experimental)','200 W (Numerical)','200 W (Num. Losses)','400 W (Experimental)','400 W (Numerical)','400 W (Num. Losses)');
% %set(h1,'FontSize',12,'Interpreter','latex');
% set(h1,'FontSize',12,'FontName','Times New Roman','Location','NorthWest');

% % for i=1:11
% %     TCout2(1,i)=TCout(1,12-i);
% % end
% % for i=1:5
% %     TCout2(3,i)=TCout(3,6-i);
% % end
% % 
% % figure1 = figure('PaperSize',[20.98 29.68]);
% % axes('Parent',figure1,'FontSize',12);
% % %xlim([280 320]);
% % %ylim([0 20]);
% % box('on');
% % grid('on');
% % hold('all');
% % plot(TCout2(1,1:11),Qloss_valve_N_dT(1,1:11),'bo','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',6)
% % plot(TCout(2,1:11),Qloss_valve_N_dT(2,1:11),'r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',6)
% % plot(TCout2(3,1:5),Qloss_valve_V_dT(1,1:5),':bo','Linewidth',1.5','MarkerSize',6)
% % plot(TCout(4,1:4),Qloss_valve_V_dT(1,1:4),':ro','Linewidth',1.5','MarkerSize',6)
% % plot(TC_F(1,:),Qloss_valve_F_dT(1,:),':co','Linewidth',1.5','MarkerSize',6)
% % % plot(TH_N(:),dTint(2,:),':r^','Linewidth',1.5','MarkerSize',6)
% % % plot(TH_N(1,1:9),dTloss_int_N(5,1:9),'--r^','Linewidth',1.5','MarkerSize',6)
% % print('-depsc','fig_DeltaT.eps');
% % xlabel('Cold end temperature (K)', 'Fontsize', 14)
% % ylabel('Valve Heat loss (W)', 'Fontsize', 14)
% % h1 = legend('200W THin dep','400W THin dep','200W V dep','400W V dep','200W F dep');
% % 
% % 
% %Volumetric flow rate dependence:
V(V==0) = nan; %plot nonzero values
%DeltaT(DeltaT==0) = nan; %plot nonzero values
%Utilization = [0.374669312169312,0.562003968253968,0.749338624338625,0.936673280423281,1.12400793650794,1.31134259259259,1.49867724867725,1.68601190476191,1.87334656084656,2.24801587301587,2.62268518518519,2.99735449735450]; %for V_V_2 = [100,150,200,250,300,350,400,450,500,600,700,800]; at 1.5 Hz (File UtilizationProtoGdLong.m) 
Utilization = [0.0936673280423281,0.140500992063492,0.187334656084656,0.234168320105820,0.281001984126984,0.327835648148148,0.374669312169312,0.421502976190476,0.468336640211640,0.562003968253968,0.655671296296296,0.749338624338625;] %for V_V_2 = [100,150,200,250,300,350,400,450,500,600,700,800]; at 1.5 Hz (File UtilizationProtoGdLong2.m) 
figure5 = figure('PaperSize',[20.98 29.68]);
%a1=axes('Parent',figure5,'XAxisLocation','top','FontSize',14);
%plot(Utilization(1,1:12),dTint_V_2(1,1:12),':b+','Linewidth',1','MarkerSize',6,'Parent',a1)
%xlabel('XAxisLocation','top','Utilization (-)', 'Fontsize', 16)
ax1=axes('Parent',figure5,'FontSize',20,'FontName','Times New Roman');
hold on
xlim([100 800]);
% set(get((a2),'XLim'),[100 800]);
% set(get((a1),'XLim'),[0.374669312169312 2.99735449735450]);
% set(get((a1),'YLim'),[0 25]);
% set(get((a2),'YLim'),[0 25]);
ylim([0 25]);
box('on');
grid('on');
hold('all');
plot(V(3,1:5),deltaT(3,1:5),'bo','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',9,'Parent',ax1)
plot(V_V_2(1:12),dTint_V_2(1,1:12),':b','Linewidth',2','MarkerSize',8,'Parent',ax1)
plot(V_V_2(1:12),dTloss_int_V_2(3,1:12),'-b','Linewidth',1.5','MarkerSize',7,'Parent',ax1)
% plot(V_V(:),dTint_V(1,:),':bo','Linewidth',1.5','MarkerSize',6)
% plot(V_V(:),dTloss_int_V(3,:),'--bo','Linewidth',1.5','MarkerSize',6)
plot(V(4,:),deltaT(4,:),'r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',9,'Parent',ax1)
plot(V_V_2(:),dTint_V_2(2,:),':r','Linewidth',2','MarkerSize',9,'Parent',ax1)
plot(V_V_2(:),dTloss_int_V_2(5,:),'-r','Linewidth',1.5','MarkerSize',8,'Parent',ax1)
% plot(V_V(:),dTint_V(2,:),':r^','Linewidth',1.5','MarkerSize',6)
% plot(V_V(:),dTloss_int_V(5,:),'--r^','Linewidth',1.5','MarkerSize',6)
print('-depsc','fig_DeltaT.eps');
xlabel('Volumetric flow rate (L/h)', 'Fontsize', 22)
%set(get((a1),'Xlabel'),'String','Utilization (-)', 'Fontsize', 16)
%set(get((a2),'Xlabel'),'String','Volumetric flow rate (L/h)', 'Fontsize', 16)
ylabel('Temperature Span (K)', 'Fontsize', 22)
h1 = legend('200 W (Experimental)','200 W (Numerical)','200 W (Num. losses)','400 W (Experimental)','400 W (Numerical)','400 W (Num. losses)');
set(h1,'FontSize',18,'Location','SouthEast','FontName','Times New Roman');

ax2 = axes('XAxisLocation','top','Color','none','XColor','k','YColor','k','Parent',figure5,'FontSize',18,'YTickLabel',[],'FontName','Times New Roman');
%'Position',get(ax1,'Position')
linkaxes([ax1 ax2],'y');
hold on
plot(Utilization(1,1:12),dTint_V_2(1,1:12),':b','Linewidth',2','MarkerSize',8,'Parent',ax2)
plot(Utilization(1,1:12),dTloss_int_V_2(3,1:12),'-b','Linewidth',1.5','MarkerSize',7,'Parent',ax2)
xlabel('Utilization factor (-)', 'Fontsize', 22,'FontName','Times New Roman')
%set(get((ax2),'XLim'));
ylim([0 25]);
xlim([0.0936673280423281 0.749338624338625]);




freq_F_2(1,1) = 0.25;
dTint_F_2(1,1) = 1.02;
dTloss_int_F_2(3,1) = 1.04;
for i=1:length(freq_F)
    freq_F_2(i+1) = freq_F(i);
    dTint_F_2(1,i+1) = dTint_F(1,i);
    dTloss_int_F_2(3,i+1) = dTloss_int_F(3,i);
end
       
fexp_400 = [0.25 0.5 1 1.5 2.25 3 4];
dTexp_400 = [3.92 11.69 15.4 16.81 17.4 16.8 15.77];
fexp_200 = [0.25 0.5 1 1.5];
dTexp_200 = [4.28 7.83 8.59 4.99];
% 
figure1 = figure('PaperSize',[20.98 29.68]);
axes('Parent',figure1,'FontSize',14);
%xlim([280 320]);
ylim([0 25]);
box('on');
grid('on');
hold('all');
plot(freq(5,:),deltaT(5,:),'ko','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',8)
%plot(freq_F_2(:),dTint_F_2(1,:),':ko','Linewidth',1.5','MarkerSize',6)
% plot(freq_F_2(:),dTloss_int_F_2(3,:),'--ko','Linewidth',1.5','MarkerSize',6)
plot(fexp_400,dTexp_400,'b^','Linewidth',2','MarkerFaceColor','b','MarkerSize',8)
plot(fexp_200,dTexp_200,'rs','Linewidth',2','MarkerFaceColor','r','MarkerSize',8)
%plot(freq_F(:),dTloss_int_F_novalve(1,:),'--co','Linewidth',1.5','MarkerSize',6)
% plot(THin(2,:),deltaT(2,:),'r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',6)
% plot(TH_N(:),dTint(2,:),':r^','Linewidth',1.5','MarkerSize',6)
% plot(TH_N(1,1:9),dTloss_int_N(5,1:9),'--r^','Linewidth',1.5','MarkerSize',6)
print('-depsc','fig_DeltaT.eps');
xlabel('Frequency (Hz)', 'Fontsize', 16)
ylabel('Temperature Span (K)', 'Fontsize', 16)
%h1 = legend('500 L/h & 200 W Experimental','500 L/h & 200 W Numerical','500 L/h & 200 W Num. losses','400 L/h - Experimental','200 L/h - Experimental');
%h1 = legend('500 L/h - Numerical','500 L/h - Num. losses','500 L/h - Experimental','400 L/h - Experimental','200 L/h - Experimental');
h1 = legend('500 L/h & 200 W Experimental','400 L/h & 200 W Experimental','200 L/h & 200 W Experimental');
set(h1,'Fontsize',13)

%% PLOT comparison between different simulations after changing the mass flow rate profile:
% 
% figure1 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure1,'FontSize',12);
% %xlim([280 320]);
% %ylim([0 20]);
% box('on');
% grid('on');
% hold('all');
% plot(THin(1,:),deltaT(1,:),'bo','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',6)
% plot(TH_N(:),dTint(1,:),':bo','Linewidth',1.5','MarkerSize',6)
% plot(TH_N_2(:),dTint_N_2(1,:),':b^','Linewidth',1.5','MarkerSize',6)
% plot(TH_N_3(:),dTint_N_3(1,:),':b*','Linewidth',1.5','MarkerSize',6)
% plot(TH_N(:),dTloss_int_N(3,:),'--bo','Linewidth',1.5','MarkerSize',6)
% plot(TH_N_2(:),dTloss_int_N_2(3,:),'--b^','Linewidth',1.5','MarkerSize',6)
% plot(TH_N_3(:),dTloss_int_N_3(3,:),'--b*','Linewidth',1.5','MarkerSize',6)
% %plot(TH_N(:),dTloss_int_N_novalve(1,:),'--co','Linewidth',1.5','MarkerSize',6)
% plot(THin(2,:),deltaT(2,:),'ro','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',6)
% plot(TH_N(:),dTint(2,:),':ro','Linewidth',1.5','MarkerSize',6)
% plot(TH_N_2(:),dTint_N_2(2,:),':r^','Linewidth',1.5','MarkerSize',6)
% plot(TH_N_3(:),dTint_N_3(2,:),':r*','Linewidth',1.5','MarkerSize',6)
% plot(TH_N(1,1:9),dTloss_int_N(5,1:9),'--ro','Linewidth',1.5','MarkerSize',6)
% plot(TH_N_2(1,1:9),dTloss_int_N_2(5,1:9),'--r^','Linewidth',1.5','MarkerSize',6)
% plot(TH_N_3(1,1:9),dTloss_int_N_3(5,1:9),'--r*','Linewidth',1.5','MarkerSize',6)
% %plot(TH_N(1,1:9),dTloss_int_N_novalve(2,1:9),'--c^','Linewidth',1.5','MarkerSize',6)
% print('-depsc','fig_DeltaT.eps');
% xlabel('Hot end in Temperature (K)', 'Fontsize', 14)
% ylabel('Temperature Span (K)', 'Fontsize', 14)
% h1 = legend('200W Exp','200W no loss old 440','200W no loss new 400','200W no loss new 440','200W loss old 440','200W loss new 400','200W loss new 440','400W Exp','400W no loss old 440','400W no loss new 400','400W no loss new 440','400W loss old 440','400W loss new 400','400W loss new 440');
% %'200W loss novalve'
% 
% %Volumetric flow rate dependence after changing flow rate profile:
% %DeltaT(DeltaT==0) = nan; %plot nonzero values
% figure5 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure5,'FontSize',12);
% %xlim([280 320]);
% %ylim([0 20]);
% box('on');
% grid('on');
% hold('all');
% plot(V(3,:),deltaT(3,:),'bo','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',6)
% plot(V_V(:),dTint_V(1,:),':bo','Linewidth',1.5','MarkerSize',6)
% plot(V_V_2(:),dTint_V_2(1,:),':b*','Linewidth',1.5','MarkerSize',6)
% plot(V_V(:),dTloss_int_V(3,:),'--bo','Linewidth',1.5','MarkerSize',6)
% plot(V_V_2(:),dTloss_int_V_2(3,:),'--b*','Linewidth',1.5','MarkerSize',6)
% plot(V(4,:),deltaT(4,:),'r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',6)
% plot(V_V(:),dTint_V(2,:),':r^','Linewidth',1.5','MarkerSize',6)
% plot(V_V_2(:),dTint_V_2(2,:),':r*','Linewidth',1.5','MarkerSize',6)
% plot(V_V(:),dTloss_int_V(5,:),'--r^','Linewidth',1.5','MarkerSize',6)
% plot(V_V_2(:),dTloss_int_V_2(5,:),'--r*','Linewidth',1.5','MarkerSize',6)
% print('-depsc','fig_DeltaT.eps');
% xlabel('Volumetric flow rate (L/h)', 'Fontsize', 14)
% ylabel('Temperature Span (K)', 'Fontsize', 14)
% h1 = legend('200W Exp','200W no loss old','200W no loss new','200W loss old','200W loss new','400W Exp','400W no loss old','400W no loss new','400W loss old','400W loss new');

%% PLOT Losses at the Hot end:

% plot(THin(1,1:11),Qloss_tube_H(1,1:11))
% plot(THin(1,1:11),Qloss_flowhead_H(1,1:11))
% plot(THin(1,1:11),Qloss_valve_H(1,1:11))
% plot(THin(1,1:11),Qloss_tot_exp_H(1,1:11))
% 
% figure1 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure1,'FontSize',12);
% %xlim([280 320]);
% %ylim([0 20]);
% box('on');
% grid('on');
% hold('all');
% plot(THin(1,1:11),Qloss_tube_H(1,1:11),'--bo','Linewidth',1.5','MarkerSize',6)
% plot(THin(1,1:11),Qloss_flowhead_H(1,1:11),'--b^','Linewidth',1.5','MarkerSize',6)
% plot(THin(1,1:11),Qloss_valve_H(1,1:11),'--b*','Linewidth',1.5','MarkerSize',6)
% plot(THin(1,1:11),Qloss_tot_exp_H(1,1:11),':ro','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',6)
% % plot(THin(2,:),deltaT(2,:),'r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',6)
% % plot(TH_N(:),dTint(2,:),':r^','Linewidth',1.5','MarkerSize',6)
% % plot(TH_N(1,1:9),dTloss_int_N(5,1:9),'--r^','Linewidth',1.5','MarkerSize',6)
% print('-depsc','fig_DeltaT.eps');
% xlabel('Hot in end temperature (K)', 'Fontsize', 14)
% ylabel('Heat loss (W)', 'Fontsize', 14)
% h1 = legend('- Q tube','- Q flowhead','Q valve','Q total');

%% PLOT Losses Bars (Kaspar):

TClosses = 275:2.5:300; % Different cold end temperatures
THlosses = 300; % Suppose that the hot end is always 5 K above room temperature
Troomlosses = 295;
freqlosses = 1.5;
%Calculate losses:
%Qloss_tubes_Cu_losses=zeros(1,length(TClosses));
%Qloss_tubes_plas_losses=zeros(1,length(TClosses));
Qloss_tubes_losses=zeros(1,length(TClosses));
Qloss_filter_losses =zeros(1,length(TClosses));
Rtubes_losses = 3; %For thermag V presentation - only the Cupper tubes were insulated. The plastic not. 
for j=1:length(TClosses)
    Qloss_tubes_losses(1,j)=(Troomlosses-TClosses(1,j))/Rtubes_Cu_av + (Troomlosses-TClosses(1,j))/Rtubes_plas_av;  %W Heat loss in the Cu tubes
    %Qloss_tubes_Cu_losses(1,j)=(Troomlosses-TClosses(1,j))/Rtubes_Cu_av;  %W Heat loss in the Cu tubes
    %Qloss_tubes_plas_losses(1,j)=(Troomlosses-TClosses(1,j))/Rtubes_plas_av;  %W Heat loss in the plastic tubes
    %Qloss_tubes_losses(1,j) = Qloss_tubes_Cu_losses(1,j) + Qloss_tubes_plas_losses(1,j);  %W Heat loss in the plastic tubes
    Qloss_filter_losses(1,j) = (Troomlosses-TClosses(1,j))/Rfilter_av; %W Heat loss in the filter
    Qloss_flowhead_losses(1,j)=(Troomlosses-TClosses(1,j))/Rflowhead; %W
end
Tvalve=10;
Qloss_valve_losses = Tvalve.*freqlosses;
%Qloss_valve_losses=13.2*freqlosses-2.5;
Qloss_reg_losses = zeros(1,length(TClosses));
for j=1:length(TClosses)
    if TClosses(1,j)<Troomlosses
            Qloss_reg_losses(1,j) = (Troomlosses-TClosses(1,j))^2 / (2*Rreg*(THlosses-TClosses(1,j)));  
        else
            Qloss_reg_losses(1,j) = (Troomlosses-TClosses(1,j))/Rreg; 
    end
end
Qloss_tot_losses = Qloss_tubes_losses + Qloss_flowhead_losses + Qloss_reg_losses + Qloss_valve_losses + Qloss_filter_losses;

%dat = importdata( 'data.txt' );
%barData = zeros( length(dat(:,1)),length(dat(1,:)) );
datalosses = zeros(5,length(TClosses(1,:)));
barData = zeros(5,length(TClosses(1,:)));

%the total bar
datalosses(1,:) = TClosses(1,:);
%datalosses(2,:) = Qloss_valve_losses(1,:);
datalosses(2,:) = Qloss_reg_losses(1,:);
datalosses(3,:) = Qloss_flowhead_losses(1,:);
%datalosses(5,:) = Qloss_tubes_Cu_losses(1,:);
datalosses(4,:) = Qloss_tubes_losses(1,:);
datalosses(5,:) = Qloss_filter_losses(1,:);

barData(1,:) = sum(datalosses(2:end,:),1);

posSum = zeros(length(TClosses(1,:)),1);
negSum = zeros(length(TClosses(1,:)),1);
%each individual bar
for i=2:length(datalosses(:,1))
	for j=1:length(datalosses(i,:))
        if datalosses(i,j) > 0             
            barData(i,j) = datalosses(i,j) + posSum(j);
            posSum(j) = posSum(j) + datalosses(i,j);
        else
            barData(i,j) = datalosses(i,j) + negSum(j);
            negSum(j) = negSum(j) + datalosses(i,j);
        end
	end
end
% 
% filename = 'lossBar';
% %names = {'Total','Valve','Filter','Tubes', 'Flowhead','Regenerator'};
% %names = {'Total','Flow distributor','Regenerator','Flowhead','Tubes','Filter'};
% names = {'Total','Regenerator','Flowhead','Tubes','Filter'};
% doGray = true;
% %doPlotBarLoss( datalosses(1,:), barData,doGray,filename, names);
% figure1 = figure('PaperSize',[20.98 29.68]);
% axes1 = axes('Parent',figure1,'YGrid','on','XGrid','on','LineWidth',2,'FontSize',14,'FontName','Times New Roman');
% xlim([270 305]);
% ylim([-25 80]);
% box(axes1,'on');
% hold(axes1,'all');
% %col = gray( length(ymat(:,1)) );
% if doGray 
%     col = gray( length(barData(:,1)) );
% else
%     col = jet( length(barData(:,1)) );
% end
% for i=length(barData(:,1)):-1:2
%     bar(datalosses(1,:),barData(i,:),'FaceColor',col(i,:),'Parent',axes1,'DisplayName',char(names(i)));
% end
% plot(datalosses(1,:),barData(1,:)+15,'ko:','Linewidth',1.5','MarkerFaceColor','w','MarkerSize',6)
% % for i=length(ymat(:,1)):-1:2
% %     bar(xvec,ymat(i,:),'FaceColor',col(i,:),'Parent',axes1,'DisplayName',char(names(i)));
% % end
% %xlabel('T_{cold} [K]','FontWeight','bold','FontSize',16);
% xlabel('Cold reservoir temperature (K)','FontSize',17,'FontName','Times New Roman');
% ylabel('Heat loss (W)','Fontsize', 17,'FontName','Times New Roman');
% %h1=legend(axes1,'show');
% %h2=legend(axes1,'Total');
% h1=legend(axes1,'$\dot{Q}_{\rm filter}$','$\dot{Q}_{\rm tubes}$','$\dot{Q}_{\rm f{}lowhead}$','$\dot{Q}_{\rm regenerator}$','$\dot{Q}_{\rm total}$');
% set(h1,'FontSize',14,'Interpreter','latex');
% %set(h2,'FontSize',14);
% %applyhatch( figure1, '/\|-+.' );
% print( '-depsc', [filename 'lossBar.eps'] );

%% PLOTS COP from experimental data
% %hot temperature dependence:
% figure2 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure2,'FontSize',19,'FontName','Times New Roman');
% xlim([285 307]);
% %ylim([0 20]);
% box('on');
% grid('on');
% hold('all');
% plot(THin(1,:),COP(1,:),'bo','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',9)
% %plot(TH_N(1,:),COPint_N(1,:),':bo','Linewidth',1.5','MarkerSize',6)
% plot(THin(1,:),COP_novalve(1,:),'--b','Linewidth',1.8','MarkerFaceColor','w','MarkerSize',9)
% %plot(THin(1,:),COP_noloss(1,:),'-b','Linewidth',1.5','MarkerSize',6)
% %plot(THin(1,:),COP_AMR(1,:),'--bo','Linewidth',1.5','MarkerSize',6)
% plot(THin(1,:),COP_motor_eff(1,:),'-b','Linewidth',1.8','MarkerFaceColor','w','MarkerSize',9)
% plot(THin(2,:),COP(2,:),'r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',9)
% %plot(TH_N(1,:),COPint_N(2,:),':r^','Linewidth',1.5','MarkerSize',6)
% plot(THin(2,:),COP_novalve(2,:),'--r','Linewidth',1.8','MarkerFaceColor','w','MarkerSize',9)
% %plot(THin(2,:),COP_noloss(2,:),'--r^','Linewidth',1.5','MarkerSize',6)
% %plot(THin(2,:),COP_AMR(2,:),'--r^','Linewidth',1.5','MarkerSize',6)
% plot(THin(2,:),COP_motor_eff(2,:),'-r','Linewidth',1.8','MarkerFaceColor','w','MarkerSize',9)
% %plot(TH_N(1,:),COP_carnot(1,:),'-r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',6)
% %plot(THin(1,:),COP_Qcalc(1,:),'-c^','Linewidth',1.5','MarkerFaceColor','c','MarkerSize',6)
% %plot(THin(1,:),COP_QCcalc(1,:),'-k^','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',6)
% print('-depsc','fig_DeltaT.eps');
% xlabel('Hot reservoir temperature (K)', 'Fontsize', 21,'FontName','Times New Roman')
% ylabel('COP (-)', 'Fontsize', 21,'FontName','Times New Roman')
% h2 = legend('200 W (Experimental)','200 W (No flow dist.)','200 W (AMR)','400 W (Experimental)','400 W (No flow dist.)','400 W (AMR)');
% %h2 = legend('200W Exp','200W no valve','200W ideal','200W Motor eff','400W Exp','400W no valve','400W ideal','400W Motor eff');
% set(h2,'FontSize',16,'FontName','Times New Roman');
% 
% % Comparison COPs:
% figure2 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure2,'FontSize',19,'FontName','Times New Roman');
% xlim([285 307]);
% %ylim([0 20]);
% box('on');
% grid('on');
% hold('all');
% plot(THin(1,:),COP_ideal(1,:),'bo','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',9)
% plot(THin(1,:),COP_carnot_Tch(1,:),'--b','Linewidth',1.8','MarkerFaceColor','w','MarkerSize',9)
% plot(THin(1,:),COP_carnot_TH(1,:),'-b','Linewidth',1.8','MarkerFaceColor','w','MarkerSize',9)
% print('-depsc','fig_DeltaT.eps');
% xlabel('Hot reservoir temperature (K)', 'Fontsize', 21,'FontName','Times New Roman')
% ylabel('COP (-)', 'Fontsize', 21,'FontName','Times New Roman')
% h2 = legend('200 W (Experimental)','200 W (No flow dist.)','200 W (AMR)','400 W (Experimental)','400 W (No flow dist.)','400 W (AMR)');
% %h2 = legend('200W Exp','200W no valve','200W ideal','200W Motor eff','400W Exp','400W no valve','400W ideal','400W Motor eff');
% set(h2,'FontSize',16,'FontName','Times New Roman');

% %hot temperature dependence (COP_cycle):
% figure2 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure2,'FontSize',19,'FontName','Times New Roman');
% xlim([285 307]);
% ylim([0 8]);
% box('on');
% grid('on');
% hold('all');
% plot(THin(1,:),COP_cycle(1,:),'-b','Linewidth',2','MarkerFaceColor','w','MarkerSize',9)
% plot(THin(2,:),COP_cycle(2,:),'-r','Linewidth',2','MarkerFaceColor','w','MarkerSize',9)
% print('-depsc','fig_DeltaT.eps');
% xlabel('Hot reservoir temperature (K)', 'Fontsize', 21,'FontName','Times New Roman')
% ylabel('COP (-)', 'Fontsize', 21,'FontName','Times New Roman')
% h2 = legend('200 W','400 W');
% %h2 = legend('200W Exp','200W no valve','200W ideal','200W Motor eff','400W Exp','400W no valve','400W ideal','400W Motor eff');
% set(h2,'FontSize',16,'Location','SouthEast','FontName','Times New Roman');

% 
% % Volumetric flow rate dependence:
% figure6 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure6,'FontSize',15,'FontName','Times New Roman');
% xlim([100 700]);
% ylim([0 3]);
% box('on');
% grid('on');
% hold('all');
% plot(V(3,:),COP(3,:),'bo','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',7)
% plot(V(3,:),COP_novalve(3,:),'--b','Linewidth',1.5','MarkerSize',6)
% %plot(V(3,:),COP_noloss(3,:),'--bo','Linewidth',1.5','MarkerSize',6)
% %plot(V_V(1,:),COPint_V(1,:),':bo','Linewidth',1.5','MarkerSize',6)
% %plot(V(3,:),COP_AMR(3,:),'--bo','Linewidth',1.5','MarkerSize',6)
% plot(V(3,:),COP_motor_eff(3,:),'-b','Linewidth',1.5','MarkerSize',6)
% plot(V(4,:),COP(4,:),'r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',6)
% plot(V(4,:),COP_novalve(4,:),'--r','Linewidth',1.5','MarkerSize',6)
% %plot(V(4,:),COP_noloss(4,:),'--r^','Linewidth',1.5','MarkerSize',6)
% %plot(V_V(1,:),COPint_V(2,:),':r^','Linewidth',1.5','MarkerSize',6)
% %plot(V(4,:),COP_AMR(4,:),'--r^','Linewidth',1.5','MarkerSize',6)
% plot(V(4,:),COP_motor_eff(4,:),'-r','Linewidth',1.5','MarkerSize',6)
% %plot(THin(1,:),COP_carnot(1,:),'-r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',6)
% %plot(THin(1,:),COP_Qcalc(1,:),'-c^','Linewidth',1.5','MarkerFaceColor','c','MarkerSize',6)
% %plot(THin(1,:),COP_QCcalc(1,:),'-k^','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',6)
% print('-depsc','fig_DeltaT.eps');
% xlabel('Volumetric flow rate (L/h)', 'Fontsize', 17,'FontName','Times New Roman')
% ylabel('COP (-)', 'Fontsize', 17,'FontName','Times New Roman')
% h2 = legend('200 W (Experimental)','200 W (No flow dist.)','200 W (AMR)','400 W (Experimental)','400 W (No flow dist.)','400 W (AMR)');
% set(h2,'FontSize',13,'FontName','Times New Roman');

fexp_400 = [0.25 0.5 1 1.5 2.25 3 4];
COPexp_400 = [1.44	1.03	0.81	0.70	0.53	0.49	0.45];
fexp_200 = [0.25 0.5 1 1.5];
COPexp_200 = [1.84	1.55	1.14	1.00];
% 
% % Frequency dependence:
% figure6 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure6,'FontSize',14);
% %xlim([150 650]);
% %ylim([0 20]);
% box('on');
% grid('on');
% hold('all');
% plot(freq(5,:),COP(5,:),'ko','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',8)
% plot(fexp_400,COPexp_400,'b^','Linewidth',2','MarkerFaceColor','b','MarkerSize',8)
% plot(fexp_200,COPexp_200,'rs','Linewidth',2','MarkerFaceColor','r','MarkerSize',8)
% %plot(V(3,:),COP_novalve(3,:),':bo','Linewidth',1.5','MarkerSize',6)
% %plot(V(3,:),COP_noloss(3,:),'--bo','Linewidth',1.5','MarkerSize',6)
% %plot(V_V(1,:),COPint_V(1,:),':bo','Linewidth',1.5','MarkerSize',6)
% %plot(V(3,:),COP_AMR(3,:),'--bo','Linewidth',1.5','MarkerSize',6)
% % plot(V(3,:),COP_motor_eff(3,:),'--bo','Linewidth',1.5','MarkerSize',6)
% % plot(V(4,:),COP(4,:),'r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',6)
% % plot(V(4,:),COP_novalve(4,:),':ro','Linewidth',1.5','MarkerSize',6)
% %plot(V(4,:),COP_noloss(4,:),'--r^','Linewidth',1.5','MarkerSize',6)
% %plot(V_V(1,:),COPint_V(2,:),':r^','Linewidth',1.5','MarkerSize',6)
% %plot(V(4,:),COP_AMR(4,:),'--r^','Linewidth',1.5','MarkerSize',6)
% % plot(V(4,:),COP_motor_eff(4,:),'--r^','Linewidth',1.5','MarkerSize',6)
% %plot(THin(1,:),COP_carnot(1,:),'-r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',6)
% %plot(THin(1,:),COP_Qcalc(1,:),'-c^','Linewidth',1.5','MarkerFaceColor','c','MarkerSize',6)
% %plot(THin(1,:),COP_QCcalc(1,:),'-k^','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',6)
% print('-depsc','fig_DeltaT.eps');
% xlabel('Frequency (Hz)', 'Fontsize', 16)
% ylabel('COP (-)', 'Fontsize', 16)
% h2 = legend('500 L/h & 200 W Experimental','400 L/h & 200 W Experimental','200 L/h & 200 W Experimental');
% set(h2,'FontSize',14);

%% PLOTS COP Carnot and Numerical:
%hot temperature dependence:

% figure7 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure7,'FontSize',12);
% %xlim([280 320]);
% %ylim([0 20]);
% box('on');
% grid('on');
% hold('all');
% plot(TH_N(1,:),COP_N(1,:),':bo','Linewidth',1.5','MarkerSize',6)
% plot(THin(1,:),COP_carnot(1,:),'bo','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',6)
% plot(TH_N(1,:),COP_N(2,:),':r^','Linewidth',1.5','MarkerSize',6)
% plot(THin(2,:),COP_carnot(2,:),'r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',6)
% %plot(THin(1,:),COP_Qcalc(1,:),'-c^','Linewidth',1.5','MarkerFaceColor','c','MarkerSize',6)
% %plot(THin(1,:),COP_QCcalc(1,:),'-k^','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',6)
% print('-depsc','fig_DeltaT.eps');
% xlabel('Hot end in Temperature (K)', 'Fontsize', 14)
% ylabel('COP (-)', 'Fontsize', 14)
% h7 = legend('COP 200W Num no loss','COP Carnot 200W Exp','COP 400W Num no loss','COP Carnot 400W Exp');
% % 
% % % % Volumetric flow rate dependence:
% COP_carnot(4,4) = NaN; %crazy COP (too low temperature span)
% 
% figure8 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure8,'FontSize',12);
% xlim([150 850]);
% %ylim([0 20]);
% box('on');
% grid('on');
% hold('all');
% plot(V_V(1,:),COP_V(1,:),':bo','Linewidth',1.5','MarkerSize',6)
% plot(V(3,:),COP_carnot(3,:),'bo','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',6)
% plot(V_V(1,:),COP_V(2,:),':r^','Linewidth',1.5','MarkerSize',6)
% plot(V(4,:),COP_carnot(4,:),'r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',6)
% %plot(THin(1,:),COP_Qcalc(1,:),'-c^','Linewidth',1.5','MarkerFaceColor','c','MarkerSize',6)
% %plot(THin(1,:),COP_QCcalc(1,:),'-k^','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',6)
% print('-depsc','fig_DeltaT.eps');
% xlabel('Volumetric flow rate (L/h)', 'Fontsize', 14)
% ylabel('COP (-)', 'Fontsize', 14)
% h8 = legend('COP 200W Num no loss','COP Carnot 200W Exp','COP 400W Num no loss','COP Carnot 400W Exp');

%% PLOTS Heat Loss:
%hot temperature dependence:

% figure3 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure3,'FontSize',10);
% %xlim([280 320]);
% %ylim([0 20]);
% box('on');
% grid('on');
% hold('all');
% plot(THin(1,:),Qloss(1,:),'-b^','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',8)
% plot(THin(2,:),Qloss(2,:),'-r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',8)
% print('-depsc','fig_DeltaT.eps');
% xlabel('Hot end in Temperature (K)', 'Fontsize', 14)
% ylabel('Heat Loss (W)', 'Fontsize', 14)
% h3 = legend('200 W','400 W');

%% PLOTS Exergetic Equivalent:

% %Hot temperature dependence:
% %with double YY:
% % figure4 = figure('PaperSize',[20.98 29.68]);
% % %axes('Parent',figure4,'FontSize',12);
% % %xlim([280 320]);
% % %ylim([0 20]);
% % box('on');
% % grid('on');
% % hold('all');
% % [AX,H1,H2]=plotyy(THin(1,:),Ex_Q(1,:),THin(1,:),eta_Ex(1,:))
% % [AX,H3,H4]=plotyy(THin(2,:),Ex_Q(2,:),THin(2,:),eta_Ex(2,:))
% % %plotyy(THin(1,:),Ex_Q(1,:),THin(1,:),eta_Ex(1,:),'plot')
% % %plotyy(THin(1,:),Ex_Q(2,:),THin(2,:),eta_Ex(2,:),'plot')
% % set(get(AX(1),'Ylabel'),'String','Exergy (W)', 'Fontsize', 14) 
% % set(get(AX(2),'Ylabel'),'String','Exergy Efficiency (-)', 'Fontsize', 14) 
% % %title('Exergy') 
% % set(H1,'LineStyle','-','Color','b','Marker','^','MarkerEdgeColor','b','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',8)
% % set(H2,'LineStyle','-','Color','b','Marker','o','MarkerEdgeColor','b','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',8)
% % set(H3,'LineStyle','-','Color','r','Marker','^','MarkerEdgeColor','r','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',8)
% % set(H4,'LineStyle','-','Color','r','Marker','o','MarkerEdgeColor','r','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',8)
% % print('-depsc','fig_DeltaT.eps');
% % xlabel('Hot end in Temperature (K)', 'Fontsize', 14)
% % %AX(1) = ylabel('Exergy (W)', 'Fontsize', 14)
% % %AX(2) = ylabel('Exergy Efficiency (-)', 'Fontsize', 14)
% % h4 = legend('200 W','400 W');
% 
% figure4 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure4,'FontSize',15,'FontName','Times New Roman');
% %xlim([280 320]);
% ylim([0 20]);
% box('on');
% grid('on');
% hold('all');
% plot(THin(1,:),Ex_Q(1,:),':bo','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',8)
% plot(THin(2,:),Ex_Q(2,:),':r^','Linewidth',1.5','MarkerFaceColor','w','MarkerSize',8)
% print('-depsc','fig_DeltaT.eps');
% xlabel('Hot reservoir temperature (K)', 'Fontsize', 15,'FontName','Times New Roman')
% ylabel('Exergetic Equivalent Cooling Power (W)', 'Fontsize', 15,'FontName','Times New Roman')
% h4 = legend('200 W','400 W');
% set(h4,'FontSize',16,'FontName','Times New Roman');
% % 
% figure4 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure4,'FontSize',15,'FontName','Times New Roman');
% xlim([100 700]);
% ylim([0 20]);
% box('on');
% grid('on');
% hold('all');
% plot(V(3,:),Ex_Q(3,:),':bo','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',8)
% plot(V(4,:),Ex_Q(4,:),':r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',8)
% print('-depsc','fig_DeltaT.eps');
% xlabel('Volumetric flow rate (L/h)', 'Fontsize', 15,'FontName','Times New Roman')
% ylabel('Exergetic Equivalent Cooling Power (W)', 'Fontsize', 15,'FontName','Times New Roman')
% h4 = legend('200 W','400 W');
% set(h4,'FontSize',16,'Location','NorthWest','FontName','Times New Roman');

%% PLOT Sensitivity (Kurt):

%Sensitivity to volume flow rate:
sens = zeros(1,5);
sens = [.75 0.875 1 1.125 1.25];
sens_V = zeros(2,length(sens));
%sens_V(1,1) = dTloss_int_V_2(3,3); % 200W, 1.5 Hz, 200 L/h, 22C Chiller (TH=297.7K)
sens_V(1,1) = dTloss_int_V_2(3,5); % 200W, 1.5 Hz, 300 L/h, 22C Chiller (TH=297.7K)
sens_V(1,2) = dTloss_int_V_2(3,6); % 200W, 1.5 Hz, 350 L/h, 22C Chiller (TH=297.7K)
sens_V(1,3) = dTloss_int_V_2(3,7); % 200W, 1.5 Hz, 400 L/h, 22C Chiller (TH=297.7K) 
sens_V(1,4) = dTloss_int_V_2(3,8); % 200W, 1.5 Hz, 450 L/h, 22C Chiller (TH=297.7K)
sens_V(1,5) = dTloss_int_V_2(3,9); % 200W, 1.5 Hz, 500 L/h, 22C Chiller (TH=297.7K)
%sens_V(1,7) = dTloss_int_V_2(3,10); % 200W, 1.5 Hz, 600 L/h, 22C Chiller (TH=297.7K)

sens_600 = [0.75 0.833 1 1.167 1.25];
%sens_V(2,1) = dTloss_int_V_2(5,7); % 400W, 1.5 Hz, 400 L/h, 22C Chiller (TH=297.7K)
sens_V(2,1) = dTloss_int_V_2(5,8); % 400W, 1.5 Hz, 450 L/h, 22C Chiller (TH=297.7K)
sens_V(2,2) = dTloss_int_V_2(5,9); % 400W, 1.5 Hz, 500 L/h, 22C Chiller (TH=297.7K)
sens_V(2,3) = dTloss_int_V_2(5,10); % 400W, 1.5 Hz, 600 L/h, 22C Chiller (TH=297.7K) 
sens_V(2,4) = dTloss_int_V_2(5,11); % 400W, 1.5 Hz, 700 L/h, 22C Chiller (TH=297.7K)
sens_V(2,5) = (dTloss_int_V_2(5,12) + dTloss_int_V_2(5,11)) / 2; % 400W, 1.5 Hz, 800 L/h, 22C Chiller (TH=297.7K)  ---> Substituir por 750 L/h
%sens_V(2,5) = dTloss_int_V_2(5,12) ;

% sens_200 = [.5 1 1.5];
% sens_V(3,2) = dTlossint_V_2(3,2); % 0 W, 1.5 Hz, 200 L/h, 22C Chiller (TH=297.7K) 
% sens_V(3,1) = dTlossint_V_2(3,1); % 0 W, 1.5 Hz, 100 L/h, 22C Chiller (TH=297.7K)
% sens_V(3,3) = dTlossint_V_2(3,3); % 0 W, 1.5 Hz, 300 L/h, 22C Chiller (TH=297.7K)

%Sensitivity to losses:
Qloss_sens_V=zeros(length(TC_V_2),length(V_V_2),length(sens)); %K
QCloss_sens_V=zeros(length(TC_V_2),length(V_V_2),length(sens)); %K
for i=1:length(TC_V_2)
    for j=1:length(V_V_2)
        for k=1:length(sens)
            Qloss_sens_V(i,j,k) = Qloss_tot_V_2(i,j) .* sens(k); 
            QCloss_sens_V(i,j,k) = QC_V_2(i,j) - Qloss_sens_V(i,j,k);
        end
    end
end
dTlossint_sens_V = zeros(4,length(sens));
%QCloss_N(1,2)=600;
for k=1:length(sens)
    dTlossint_sens_V(1,k) = interp1(QCloss_sens_V(:,7,k),dT_V_2(:),200); %Temperature span at 400 L/h, 1.5Hz, 22C Chiller
    dTlossint_sens_V(2,k) = interp1(QCloss_sens_V(:,7,k),dT_V_2(:),400); %Temperature span at 400 L/h, 1.5Hz, 22C Chiller
end
for k=1:length(sens)
    dTlossint_sens_V(3,k) = interp1(QCloss_sens_V(:,10,k),dT_V_2(:),200); %Temperature span at 600 L/h, 1.5Hz, 22C Chiller
    dTlossint_sens_V(4,k) = interp1(QCloss_sens_V(:,10,k),dT_V_2(:),400); %Temperature span at 600 L/h, 1.5Hz, 22C Chiller
end
 
%Sensitivity for PARTICLE SIZE and losses
sens_dp_x = [.75 .833 1 1.167 1.25];
sens_dp = zeros(2,length(sens_dp_x));
sens_dp(1,1) = dTloss_int_V_045(3,7); % 200W, 1.5 Hz, 400 L/h, 0.45 mm, 22C Chiller (TH=297.7K)
sens_dp(1,2) = dTloss_int_V_05(3,7); % 200W, 1.5 Hz, 400 L/h, 0.5 mm, 22C Chiller (TH=297.7K)
sens_dp(1,3) = dTloss_int_V_2(3,7); % 200W, 1.5 Hz, 400 L/h, 0.6 mm, 22C Chiller (TH=297.7K) 
sens_dp(1,4) = dTloss_int_V_07(3,7); % 200W, 1.5 Hz, 400 L/h, 0.7 mm, 22C Chiller (TH=297.7K)
sens_dp(1,5) = dTloss_int_V_075(3,7); % 200W, 1.5 Hz, 400 L/h, 0.75 mm, 22C Chiller (TH=297.7K)

sens_dp(2,1) = dTloss_int_V_045(5,10); % 400W, 1.5 Hz, 600 L/h, 0.45 mm, 22C Chiller (TH=297.7K)
sens_dp(2,2) = dTloss_int_V_05(5,10); % 400W, 1.5 Hz, 600 L/h, 0.5 mm, 22C Chiller (TH=297.7K)
sens_dp(2,3) = dTloss_int_V_2(5,10); % 400W, 1.5 Hz, 600 L/h, 0.6 mm, 22C Chiller (TH=297.7K) 
sens_dp(2,4) = dTloss_int_V_07(5,10); % 400W, 1.5 Hz, 600 L/h, 0.7 mm, 22C Chiller (TH=297.7K)
sens_dp(2,5) = dTloss_int_V_075(5,10); % 400W, 1.5 Hz, 600 L/h, 0.75 mm, 22C Chiller (TH=297.7K)

% figure10 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure10,'FontSize',16,'FontName','Times New Roman');
% xlim([.7 1.3]);
% %ylim([0 23]);
% box('on');
% grid('on');
% hold('all');
% %For 0 W and 200L/h:
% %plot(sens_200(1,:),sens_V(3,1:3),':co','Linewidth',1.5','MarkerFaceColor','c','MarkerSize',6)
% %plot(sens(1,:),dTlossint_sens_V(1,:),':b^','Linewidth',1.5','MarkerSize',6)
% 
% %For 200 W and 400 L/h:
% plot(sens_dp_x(1,:),sens_dp(1,:),':bs','Linewidth',1.4','MarkerFaceColor','b','MarkerSize',7)
% plot(sens(1,:),dTlossint_sens_V(1,:),':b^','Linewidth',1.4','MarkerFaceColor','b','MarkerSize',7)
% plot(sens(1,:),sens_V(1,:),':bo','Linewidth',1.4','MarkerFaceColor','b','MarkerSize',7)
% 
% %For 400 W and 600L/h:
% plot(sens_dp_x(1,:),sens_dp(2,:),':rs','Linewidth',1.4','MarkerFaceColor','w','MarkerSize',7)
% plot(sens(1,:),dTlossint_sens_V(4,:),':r^','Linewidth',1.4','MarkerFaceColor','w','MarkerSize',7)
% plot(sens_600(1,:),sens_V(2,:),':ro','Linewidth',1.4','MarkerFaceColor','w','MarkerSize',7)
% 
% print('-depsc','fig_DeltaT.eps');
% xlabel('Scaling Factor (-)', 'Fontsize', 18,'FontName','Times New Roman')
% ylabel('Temperature Span (K)', 'Fontsize', 18,'FontName','Times New Roman')
% h4 = legend('200 W - 400 L/h - Sphere diameter','200 W - 400 L/h - Losses','200 W - 400 L/h - Vol. flow rate','400 W - 600 L/h - Sphere diameter','400 W - 600 L/h - Losses','400 W - 600 L/h - Vol. flow rate');
% set(h4,'FontSize',14,'Location','SouthEast','FontName','Times New Roman');


%% PLOT motor power:

MotorPower = zeros(4,5);
xMat = [0.5 1 1.5 2];
MotorPower(1,:) = [107 88 80 19 0];
MotorPower(2,:) = [160 120 105 40 0];
MotorPower(3,:) = [220 154 133 66 0];
MotorPower(4,:) = [290 192 163 127 0];
MotorPower(1:4,5) = MotorPower(1:4,2) - MotorPower(1:4,3);

% filename = 'MotorPower';
% names = {'$\dot{W}_{\rm M}$','$\dot{W}_{\rm no-f{}l}$','$\dot{W}_{\rm nl}$','$\dot{W}_{\rm f{}l} / \eta_{\rm M}$','$\dot{W}_{\rm mag} / \eta_{\rm M}$'};
% doGray = true;
% figure11 = figure('PaperSize',[20.98 29.68]);
% axes1 = axes('Parent',figure11,'YGrid','on','XGrid','on','LineWidth',2,'FontSize',16,'FontName','Times New Roman');
% %xlim([270 305]);
% box(axes1,'on');
% hold(axes1,'all');
% yMat = MotorPower(1:4,:);
% bar1 = bar(xMat,yMat,'Parent',axes1);
% if doGray 
%     col = gray( 6 );
% else
%     col = jet( 6 );
% end
% for i=1:length(bar1)
%     set(bar1(i),'DisplayName',char(names(i)),'FaceColor',col(i+1,:));
% end
% colormap gray
% 
% xlabel('Frequency (Hz)','FontSize',18,'FontName','Times New Roman');
% ylabel('Motor power (W)','Fontsize', 18,'FontName','Times New Roman');
% 
% h1=legend(axes1,'show');
% set(h1,'FontSize',14,'Interpreter','latex','Location','NorthWest');
% print('-depsc', [filename 'MotorPower.eps'] );

%% PLOTS 2nd law efficiency 
%hot temperature dependence:
figure2 = figure('PaperSize',[20.98 29.68]);
axes('Parent',figure2,'FontSize',19,'FontName','Times New Roman');
xlim([285 307]);
ylim([0 35]);
box('on');
grid('on');
hold('all');
plot(THin(1,:),eta_system(1,:)*100,'bo','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',9)
%plot(TH_N(1,:),COPint_N(1,:),':bo','Linewidth',1.5','MarkerSize',6)
plot(THin(1,:),eta_cycle(1,:)*100,':bo','Linewidth',1.8','MarkerFaceColor','w','MarkerSize',9)
%plot(THin(1,:),COP_noloss(1,:),'-b','Linewidth',1.5','MarkerSize',6)
%plot(THin(1,:),COP_AMR(1,:),'--bo','Linewidth',1.5','MarkerSize',6)
%plot(THin(1,:),COP_motor_eff(1,:),'-b','Linewidth',1.8','MarkerFaceColor','w','MarkerSize',9)
plot(THin(2,:),eta_system(2,:)*100,'r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',9)
%plot(TH_N(1,:),COPint_N(2,:),':r^','Linewidth',1.5','MarkerSize',6)
plot(THin(2,:),eta_cycle(2,:)*100,':r^','Linewidth',1.8','MarkerFaceColor','w','MarkerSize',9)
%plot(THin(2,:),COP_noloss(2,:),'--r^','Linewidth',1.5','MarkerSize',6)
%plot(THin(2,:),COP_AMR(2,:),'--r^','Linewidth',1.5','MarkerSize',6)
%plot(THin(2,:),COP_motor_eff(2,:),'-r','Linewidth',1.8','MarkerFaceColor','w','MarkerSize',9)
%plot(TH_N(1,:),COP_carnot(1,:),'-r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',6)
%plot(THin(1,:),COP_Qcalc(1,:),'-c^','Linewidth',1.5','MarkerFaceColor','c','MarkerSize',6)
%plot(THin(1,:),COP_QCcalc(1,:),'-k^','Linewidth',1.5','MarkerFaceColor','k','MarkerSize',6)
print('-depsc','fig_DeltaT.eps');
xlabel('Hot reservoir temperature (K)', 'Fontsize', 21,'FontName','Times New Roman')
ylabel('Second law efficiency (%)', 'Fontsize', 21,'FontName','Times New Roman')
h2 = legend('$\eta_{\rm 2nd}$~for $\dot{Q}_{\rm C} =$ 200 W','$\eta_{\rm 2nd,cy}$~for $\dot{Q}_{\rm C} =$ 200 W','$\eta_{\rm 2nd}$~for $\dot{Q}_{\rm C} =$ 400 W','$\eta_{\rm 2nd,cy}$~for $\dot{Q}_{\rm C} =$ 400 W');
set(h2,'FontSize',15,'Interpreter','latex');
% % 
% %volumetric flow rate dependence:
% figure2 = figure('PaperSize',[20.98 29.68]);
% axes('Parent',figure2,'FontSize',19,'FontName','Times New Roman');
% xlim([100 700]);
% ylim([0 35]);
% box('on');
% grid('on');
% hold('all');
% plot(V(3,:),eta_system(3,:)*100,'bo','Linewidth',1.5','MarkerFaceColor','b','MarkerSize',9)
% plot(V(3,:),eta_cycle(3,:)*100,':bo','Linewidth',1.8','MarkerFaceColor','w','MarkerSize',9)
% plot(V(4,:),eta_system(4,:)*100,'r^','Linewidth',1.5','MarkerFaceColor','r','MarkerSize',9)
% plot(V(4,:),eta_cycle(4,:)*100,':r^','Linewidth',1.8','MarkerFaceColor','w','MarkerSize',9)
% print('-depsc','fig_DeltaT.eps');
% xlabel('Volumetric flow rate (L/h)', 'Fontsize', 21,'FontName','Times New Roman')
% ylabel('Second law efficiency (%)', 'Fontsize', 21,'FontName','Times New Roman')
% %h2 = legend('$\eta_{system}$ for $\dot{\mathrm{Q}}_{C} =$ 200 W','$\eta_{thermod}$ for $\dot{\mathrm{Q}}_{C} =$ 200 W','$\eta_{system}$ for $\dot{\mathrm{Q}}_{C} =$ 400 W','$\eta_{thermod}$ for $\dot{\mathrm{Q}}_{C} =$ 400 W');
% h2 = legend('$\eta_{\rm 2nd}$~for $\dot{Q}_{\rm C} =$ 200 W','$\eta_{\rm 2nd,cy}$~for $\dot{Q}_{\rm C} =$ 200 W','$\eta_{\rm 2nd}$~for $\dot{Q}_{\rm C} =$ 400 W','$\eta_{\rm 2nd,cy}$~for $\dot{Q}_{\rm C} =$ 400 W');
% set(h2,'FontSize',15,'Interpreter','latex');
