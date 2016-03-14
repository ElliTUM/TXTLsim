function sol= main_Komp_2(I,I2)
% Daten laden
load('testData7')
data=data-145;
data=smooth(data,5)';
for i=1:30
data(:,i)=[];
end
global bestparam besterr
% c_0 Anfangsparameter festlegen
c_0(1)=5; % TS1-DNA in nM
c_0(2)=0; % DRp
c_0(3)=30; %RNA-Polymerase Rp in nM
c_0(4)=0; % TS1 in nM
c_0(5)=10; % Trigger DNA DTrp in nM
c_0(6)=0; % Trigger in nM
c_0(7)=0; % DTrpRp in nM
c_0(8)=10; % Anti-Trigger-DNA DTrm
c_0(9)=0; % DTrmRp
c_0(10)=0; % Anti-Trigger
c_0(11)=200; % Ribosom
c_0(12)=0; % TS1R
c_0(13)=0; % Punfold
c_0(14)=0; % P also mCer
c_0(15)=0; % TS1Tr
c_0(16)=0; % TS1TrR
c_0(17)=10; % Xm
c_0(18)=0; % TS1Xm
c_0(19)=0; % TrXm in nM
c_0(20)=0; % AntiTrXm in nM
c_0(21)=0; % TrAntiTr


% Lac-Repressor
R=2490; %Lac-Repressor-Konzentration in nM in DH5alpha Z1
%I=1*10^6; % IPTG-Konzentration in nM

KAR=0.1; %Dissoziationsrate Repressor/DNA-Bindestelle in nM
KRI=1000; %Dissoziationsrate Repressor/Inducer in nM
n1=2; % Hill-Koeffizient 
% Tet-Repressor
R2=5810; % tet-Repressor-Konzentration in nM in DH5alpha Z1
%I2=216*10^3; % aTc-Konzentration 100ng/ml

KAR2=0.1; %Dissoziationsrate Repressor/DNA-Bindestelle
KRI2=1000; %Dissoziationsrate Repressor/Inducer
n2=2; % Hill-Koeffizient

% IPTG und aTc modifizieren die Transkriptionsraten k2b und k2c

k2=1/73;
k2b=k2*1/(1+(R/KAR*(1/(1+(I/KRI)^n1))));

k2c=k2*1/(1+(R2/KAR2*(1/(1+(I2/KRI2)^n2))));


% Startparameter festlegen
param(1)=3*10^-2; % RNAP->DNA /nM /s
param(2)=5; % RD-> R + D \s Source: Belintsev et al., NRA, 1980
param(3)=1/800; % kcat /s Source: Karzbrun et al. PRL, 2011
param(4)=3*10^-2; % RNAP->DNA (Trigger-DNA) /nM /s
% param(5)=1/900; % kcat /s Source: Karzbrun et al. PRL, 2011
param(5)=k2b;
param(6)=3*10^-2; % RNAP->DNA (Anti-Trigger-DNA) /nM /s
% param(7)=1/900; % kcat /s Source: Karzbrun et al. PRL, 2011
param(7)=k2c;
param(8)=(1/100)*2.6*10^-3; % /nM/s  leaky binding rate of R
param(9)=0.015; %diss rate of R
param(10)=5/883; %translation rate
param(11)=1/90; %kfold
param(12)=2.6*10^-3; % /nM/s  binding rate of R
param(13)=7*10^-4; % /nM/s Xm association
param(14)=0.125; % /s Xm dissociation
param(15)=1/(12*60); % degradation rate
param(16)=10^-3; %Hybridisierungsrate /nM/s
param(17)=10^-5; %Dissoziationsrate TS1Tr1 geraten
param(18)=10^-3; %Hybridisierungsrate /nM/s
param(19)=10^-5; %Dissoziationsrate TS1Tr1 geraten
param(20)=10^-3; %Hybridisierungsrate /nM/s
param(21)=10^-5; %Dissoziationsrate TS1Tr1 geraten

errW=[1,1]; %error Weights [lsqE, fftE]

cutFFT=7; 
%besterr festlegen
%totalError=0; %reset Error
besterr=Inf;
bestparam=param;

% ode23s_solver aufrufen
[totalError,simData14]=ode23s_solver_Komp_2(param,c_0,data,errW,cutFFT)
sol=[simData14]
% Fitfuntion aufrufen
% Fitfkt_BS(param,c_0,data,errW,cutFFT)

