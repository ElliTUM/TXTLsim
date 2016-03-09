function sol=define_Komp(t,c,k)

%parameters
%% rate constants
kf1=k(1); kr1=k(2); k2=k(3); kf1b=k(4); k2b=k(5); kf1c=k(6); k2c=k(7); kf3=k(8); kr3=k(9); k4=k(10); kfold=k(11);
kf3b=k(12); kf5=k(13); kr5=k(14); kdeg=k(15); kf10=k(16); kr10=k(17); kf11=k(18);kr11=k(19); kf12=k(20); kr12=k(21); 


%% concentrations
%R13
D=c(1); DRp=c(2); Rp=c(3); TS1=c(4); DTrp=c(5); Tr=c(6); DTrpRp=c(7); DTrm=c(8); DTrmRp=c(9); AntiTr=c(10);

R=c(11); TS1R=c(12); Punf=c(13); P=c(14); TS1Tr=c(15); TS1TrR=c(16); 

Xm=c(17); TS1Xm=c(18); TrXm=c(19); AntiTrXm=c(20); TrAntiTr=c(21);

%% ordinary differential equations T3R5

%Transkription
%Toholdswitch Transkription, Trigger-Transkription, Anti-Trigger-Transkription
Ddot=-kf1*D*Rp+kr1*DRp+k2*DRp;
DTrpdot=-kf1b*DTrp*Rp+kr1*DTrpRp+k2b*DTrpRp;
DTrmdot=-kf1c*DTrm*Rp+kr1*DTrmRp+k2c*DTrmRp;
DRpdot=-Ddot;
DTrpRpdot=-DTrpdot;
DTrmRpdot=-DTrmdot;

Rpdot=Ddot+DTrpdot+DTrmdot;


% Tr1dot=k7*D2P-Tr1Xmdot-kdeg*Tr1Xm-TS1Tr1dot;
%leakyTranslation=Translation mit neuer Rate für kf3

%Translation und Proteinfaltung 

TS1Rdot=kf3*TS1*R-kr3*TS1R-k4*TS1R;
TS1TrRdot=kf3b*TS1Tr*R-kr3*TS1TrR-k4*TS1TrR;
Rdot=-TS1Rdot-TS1TrRdot;

Punfdot=k4*TS1R-kfold*Punf+k4*TS1TrR;
Pdot=kfold*Punf;
% 
% %Degradation für TS1,Tr1,TS2,Tr2
Xmdot=-kf5*TS1*Xm+kr5*TS1Xm+kdeg*TS1Xm-kf5*Tr*Xm+kr5*TrXm+kdeg*TrXm-kf5*AntiTr*Xm+kr5*AntiTrXm+kdeg*AntiTrXm;
TS1Xmdot=kf5*TS1*Xm-kr5*TS1Xm-kdeg*TS1Xm;
TrXmdot=kf5*Tr*Xm-kr5*TrXm-kdeg*TrXm;
AntiTrXmdot=-Xmdot-TS1Xmdot-TrXmdot;
% 


% %Zusammenlagerung TS1 und Tr1
TrAntiTrdot=kf10*Tr*AntiTr-kr10*TrAntiTr+kf11*TS1Tr*AntiTr-kr11*TS1*TrAntiTr;
TS1Trdot=kf12*TS1*Tr-kr12*TS1Tr-kf3b*TS1Tr*R+kr3*TS1TrR;
% TS1Tr1dot=kf10*TS1*Tr1-kr10*TS1*Tr1-k8*TS1Tr1*R+kr3*TS1Tr1R-kf12*TS1Tr1*TS2+kr12*TS2Tr1*TS1;
% %Zusammenlagerung TS2 und Tr2
% TS2Tr2dot=kf10*TS2*Tr2-kr10*TS2*Tr2-k8b*TS2Tr2*R+kr3*TS2Tr2R-kf12b*TS2Tr2*TS1+kr12b*TS1Tr2*TS2;
% 
% % TS1-mRNA bindet freien Trigger 2
% % TS1-mRNA ist zum gesamten Trigger komplementaer -> strand displacement vom Trigger-TS-Komplex
% % Zusammenlagerung TS1 und Tr2
% TS1Tr2dot=kf11*TS1*Tr2-kr11*TS1Tr2+kf12b*TS2Tr2*TS1-kr12b*TS1Tr2*TS2;
% %Zusammenlagerung TS2 und Tr1
% TS2Tr1dot=kf11b*TS2*Tr2-kr11b*TS2Tr1+kf12*TS1Tr1*TS2-kr12*TS2Tr1*TS1;
% 
% TS1dot=k2*DRp-TS1Rdot-TS1Xmdot-kdeg*TS1Xm-TS1Tr1-k8*TS1Tr1*R+kr3*TS1Tr1R-kf12b*TS2Tr2*TS1+kr12b*TS1Tr2*TS2;
% TS2dot=k2b*DbRp-TS2Rdot-TS2Xmdot-kdeg*TS2Xm--TS2Tr2-k8b*TS2Tr2*R+kr3*TS2Tr2R-kf12*TS1Tr1*TS2+kr12*TS2Tr1*TS1;
% 
% %Trigger-Transkription durch T7-Polymerase=P 
% D2dot=-kf6*D2*P+kr6*D2P+k2*D2P;
% Pdot=kfold*Punf+D2dot;
% D2Pdot=-D2dot;
% Tr1dot=k7*D2P-Tr1Xmdot-kdeg*Tr1Xm-TS1Tr1dot;
% 
% %% ordinary differential equations R13
% 
% %Trigger-Transkription durch T7-R13-Polymerase=Pb 
% D2bdot=-kf6b*D2b*Pb+kr6*D2bPb+k2b*D2bPb;
% Pbdot=kfold*Punfb+D2bdot;
% D2bPbdot=-D2bdot;
% Tr2dot=k7*D2bPb-Tr2Xmdot-kdeg*Tr2Xm-TS2Tr2dot;
% 
TS1dot=k2*DRp-TS1Xmdot-kdeg*TS1Xm+kf11*TS1Tr*AntiTr-kr11*TS1*TrAntiTr-kf12*TS1*Tr+kr12*TS1Tr;
Trdot=k2b*DTrpRp-TrXmdot-kdeg*TrXm-kf10*Tr*AntiTr+kr10*TrAntiTr-kf12*TS1*Tr+kr12*TS1Tr;
AntiTrdot=k2c*DTrmRp-TrAntiTrdot-AntiTrXmdot-kdeg*AntiTrXm;

sol=[Ddot;DRpdot;Rpdot;TS1dot;DTrpdot;Trdot;DTrpRpdot;DTrmdot;DTrmRpdot;AntiTrdot;Rdot;TS1Rdot;Punfdot;Pdot;TS1Trdot;TS1TrRdot;Xmdot;TS1Xmdot;TrXmdot;AntiTrXmdot;TrAntiTrdot];