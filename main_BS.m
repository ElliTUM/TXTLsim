
% Daten laden
load('testData7')
data=data-145;
data=smooth(data,5)';
for i=1:30
data(:,i)=[];
end
global bestparam besterr
% c_0 Anfangsparameter festlegen
c_0(1)=10; % DNA in nM
c_0(2)=0; %DRp
c_0(3)=30; %Rp in nM
c_0(4)=0; % mRNA in nM
c_0(5)=200; % R in nM
c_0(6)=0; % mRNAR in nM
c_0(7)=0; % Punfold 
c_0(8)=0; % P
c_0(9)=4; % Xm
c_0(10)=0; % mRNAXm
c_0(11)=5; % DNA 2 in nM
c_0(12)=0; %D2P
c_0(13)=0; % m2 in nM
c_0(14)=0; % m2Xm
c_0(15)=0; % mRNAm2
c_0(16)=0; % mRNAm2R

c_0(17)=5; % DNA in nM
c_0(18)=0; % DRp
c_0(19)=0; % mRNA in nM
c_0(20)=0; % mRNAR in nM
c_0(21)=0; % Punfold 
c_0(22)=0; % P
c_0(23)=0; % mRNAXm
c_0(24)=3; % DNA 2 in nM
c_0(25)=0; % D2P
c_0(26)=0; % m2 in nM
c_0(27)=0; % m2Xm
c_0(28)=0; % mRNAm2
c_0(29)=0; % mRNAm2R
c_0(30)=0; % mRNAm2b
c_0(31)=0; % mRNAbm2
% Startparameter festlegen

param(1)=3*10^-2; % RNAP->DNA /nM /s

param(2)=5; % RD-> R + D \s Source: Belintsev et al., NRA, 1980

param(3)=1/900; % kcat /s Source: Karzbrun et al. PRL, 2011

param(4)=(1/200)*2.6*10^-3; % /nM/s  leaky binding rate of R

param(5)=0.015; %diss rate of R

param(6)=5/883; %translation rate

param(7)=1/90; %kfold

param(8)=7*10^-4; % /nM/s Xm association

param(9)=0.125; % /s Xm dissociation

param(10)=1/(12*60); % degradation rate

param(11)=1.9*10^-4; % /nM/s T7 Polymerase ass. rate 

param(12)=0.0446; % /s T7 polymerase diss.rate

param(13)=0.0186; % kcat T7 polymerase

param(14)=2.6*10^-3; % /nM/s   binding rate of R

param(15)=10^-3; %Hybridisierungsrate /nM/s lt. Simmel-Vorlesung

param(16)=10^-5; %Dissoziationsrate TS1Tr1 geraten

param(17)=1/500; % kcat /s Source: Karzbrun et al. PRL, 2011

param(18)=(1/100)*2.6*10^-3; % /nM/s  leaky binding rate of R

param(19)=1.9*10^-4; % /nM/s T7 Polymerase ass. rate 

param(20)=1.2*10^-3; % /nM/s   binding rate of R

param(21)=1; %Hybridisierungsrate /nM/s 

param(22)=10^-5; %Dissoziationsrate TS1Tr1 geraten

param(23)=10^-8; %Hybridisierungsrate /nM/s

param(24)=10^-5; %Dissoziationsrate TS1Tr1 geraten

param(25)=10^-1; %Hybridisierungsrate /nM/s 

param(26)=10^-5; %Dissoziationsrate TS1Tr1 geraten

param(27)=10^-6; %Hybridisierungsrate /nM/s

param(28)=10^-4; %Dissoziationsrate TS1Tr1 geraten

param(29)=10^-1; %Hybridisierungsrate /nM/s 

param(30)=10^-5; %Dissoziationsrate TS1Tr1 geraten


errW=[1,1]; %error Weights [lsqE, fftE]

cutFFT=7; 
%besterr festlegen
%totalError=0; %reset Error
besterr=Inf;
bestparam=param;

% ode23s_solver aufrufen
ode23s_solver_BS(param,c_0,data,errW,cutFFT)

% Fitfuntion aufrufen
% Fitfkt_BS(param,c_0,data,errW,cutFFT)
