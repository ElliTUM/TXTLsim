function sol=define_BS(t,c,k)

%parameters
%% rate constants
kf1=k(1); kr1=k(2); k2=k(3); kf3=k(4); kr3=k(5); k4=k(6); kfold=k(7); kf5=k(8); kr5=k(9); kdeg=k(10);
kf6=k(11); kr6=k(12); k7=k(13); k8=k(14); kf10=k(15); kr10=k(16); 

k2b=k(17); kf3b=k(18);kf6b=k(19); k8b=k(20); kf10b=k(21); kr10b=k(22);kf11=k(23); kr11=k(24); kf11b=k(25);
kr11b=k(26);kf12=k(27); kr12=k(28); kf12b=k(29); kr12b=k(30);

%% concentrations
%R13
D=c(1); DRp=c(2); Rp=c(3); TS1=c(4); R=c(5); TS1R=c(6); Punf=c(7); P=c(8); Xm=c(9);TS1Xm=c(10);
D2=c(11); D2P=c(12); Tr1=c(13); Tr1Xm=c(14); TS1Tr1=c(15); TS1Tr1R=c(16);

%T3R5
Db=c(17); DbRp=c(18); TS2=c(19); TS2R=c(20); Punfb=c(21); Pb=c(22); TS2Xm=c(23);
D2b=c(24); D2bPb=c(25); Tr2=c(26); Tr2Xm=c(27); TS2Tr2=c(28); TS2Tr2R=c(29); TS1Tr2=c(30); TS2Tr1=c(31);

%% ordinary differential equations T3R5

%Transkription
Ddot=-kf1*D*Rp+kr1*DRp+k2*DRp;
Dbdot=-kf1*Db*Rp+kr1*DbRp+k2b*DbRp;
DRpdot=-Ddot;
DbRpdot=-Dbdot;

Rpdot=Ddot+Dbdot;

%leakyTranslation=Translation mit neuer Rate für kf3

%Translation und Proteinfaltung 
Rdot=-kf3*TS1*R+kr3*TS1R+k4*TS1R-k8*TS1Tr1*R+kr3*TS1Tr1R+k4*TS1Tr1R-kf3b*TS2*R+kr3*TS2R+k4*TS2R-k8b*TS2Tr2*R+kr3*TS2Tr2R+k4*TS2Tr2R;
TS1Rdot=kf3*TS1*R-kr3*TS1R-k4*TS1R;
TS2Rdot=kf3b*TS2*R-kr3*TS2R-k4*TS2R;
TS1Tr1Rdot=k8*TS1Tr1*R-kr3*TS1Tr1R-k4*TS1Tr1R;
TS2Tr2Rdot=-Rdot-TS2Rdot-TS1Tr1Rdot-TS1Rdot;

Punfdot=k4*TS1R-kfold*Punf+k4*TS1Tr1R;

Punfbdot=k4*TS2R-kfold*Punfb+k4*TS2Tr2R;

%Degradation für TS1,Tr1,TS2,Tr2
Xmdot=-kf5*TS1*Xm+kr5*TS1Xm+kdeg*TS1Xm-kf5*Tr1*Xm+kr5*Tr1Xm+kdeg*Tr1Xm-kf5*TS2*Xm+kr5*TS2Xm+kdeg*TS2Xm-kf5*Tr2*Xm+kr5*Tr2Xm+kdeg*Tr2Xm;
TS1Xmdot=kf5*TS1*Xm-kr5*TS1Xm-kdeg*TS1Xm;
TS2Xmdot=kf5*TS2*Xm-kr5*TS2Xm-kdeg*TS2Xm;
Tr1Xmdot=kf5*Tr1*Xm-kr5*Tr1Xm-kdeg*Tr1Xm;
Tr2Xmdot=-Xmdot-TS2Xmdot-TS1Xmdot-Tr1Xmdot;

%Zusammenlagerung TS1 und Tr1
TS1Tr1dot=kf10*TS1*Tr1-kr10*TS1*Tr1-k8*TS1Tr1*R+kr3*TS1Tr1R-kf12*TS1Tr1*TS2+kr12*TS2Tr1*TS1;
%Zusammenlagerung TS2 und Tr2
TS2Tr2dot=kf10*TS2*Tr2-kr10*TS2*Tr2-k8b*TS2Tr2*R+kr3*TS2Tr2R-kf12b*TS2Tr2*TS1+kr12b*TS1Tr2*TS2;

% TS1-mRNA bindet freien Trigger 2
% TS1-mRNA ist zum gesamten Trigger komplementaer -> strand displacement vom Trigger-TS-Komplex
% Zusammenlagerung TS1 und Tr2
TS1Tr2dot=kf11*TS1*Tr2-kr11*TS1Tr2+kf12b*TS2Tr2*TS1-kr12b*TS1Tr2*TS2;
%Zusammenlagerung TS2 und Tr1
TS2Tr1dot=kf11b*TS2*Tr2-kr11b*TS2Tr1+kf12*TS1Tr1*TS2-kr12*TS2Tr1*TS1;

TS1dot=k2*DRp-TS1Rdot-TS1Xmdot-kdeg*TS1Xm-TS1Tr1-k8*TS1Tr1*R+kr3*TS1Tr1R-kf12b*TS2Tr2*TS1+kr12b*TS1Tr2*TS2;
TS2dot=k2b*DbRp-TS2Rdot-TS2Xmdot-kdeg*TS2Xm--TS2Tr2-k8b*TS2Tr2*R+kr3*TS2Tr2R-kf12*TS1Tr1*TS2+kr12*TS2Tr1*TS1;

%Trigger-Transkription durch T7-Polymerase=P 
D2dot=-kf6*D2*P+kr6*D2P+k2*D2P;
Pdot=kfold*Punf+D2dot;
D2Pdot=-D2dot;
Tr1dot=k7*D2P-Tr1Xmdot-kdeg*Tr1Xm-TS1Tr1dot;

%% ordinary differential equations R13

%Trigger-Transkription durch T7-R13-Polymerase=Pb 
D2bdot=-kf6b*D2b*Pb+kr6*D2bPb+k2b*D2bPb;
Pbdot=kfold*Punfb+D2bdot;
D2bPbdot=-D2bdot;
Tr2dot=k7*D2bPb-Tr2Xmdot-kdeg*Tr2Xm-TS2Tr2dot;



sol=[Ddot;DRpdot;Rpdot;TS1dot;Rdot;TS1Rdot;Punfdot;Pdot;Xmdot;TS1Xmdot;D2dot;D2Pdot;Tr1dot;Tr1Xmdot;TS1Tr1dot;TS1Tr1Rdot;Dbdot;DbRpdot;TS2dot;TS2Rdot;Punfbdot;Pbdot;TS2Xmdot;D2bdot;D2bPbdot;Tr2dot;Tr2Xmdot;TS2Tr2dot;TS2Tr2Rdot;TS1Tr2dot;TS2Tr1dot];