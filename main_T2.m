
% Daten laden
load('testData7')
data=data-140;
data=smooth(data,5)';
for i=1:30
data(:,i)=[];
end
global bestparam besterr
% c_0 Anfangsparameter festlegen
c_0(1)=5; % DNA in nM
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
c_0(15)=0; % m2R
c_0(16)=0; % P2unf
c_0(17)=0; % P2
% Startparameter festlegen

param(1)=3*10^-2; % RNAP->DNA /nM /s

param(2)=5; % RD-> R + D \s Source: Belintsev et al., NRA, 1980

param(3)=1/900; % kcat /s Source: Karzbrun et al. PRL, 2011

param(4)=(1/20)*2.6*10^-3; % /nM/s  leaky binding rate of R

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

param(15)=5/240; % translation rate

param(16)=1/90; %kfold2

% param=[0.283569349425548 0.519568747824917 0.00331941503709697 0.000825192652206199 0.0344859540915462 0.0139527613302360 0.0265778426625039 0.00699999965010302 0.0128683835681360 0.0138850596569748 0.00166151616089664 0.00446000006871836 0.00186155191026500 0.000349890142322313 0.0105159305624588 0.106476930973575];
errW=[1,10]; %error Weights [lsqE, fftE]

cutFFT=7; 
%besterr festlegen
%totalError=0; %reset Error
besterr=Inf;
bestparam=param;

% ode23s_solver aufrufen
ode23s_solver_T2(param,c_0,data,errW,cutFFT)

% Fitfuntion aufrufen
Fitfkt_T2(param,c_0,data,errW,cutFFT)


