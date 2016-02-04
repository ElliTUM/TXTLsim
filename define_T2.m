function sol=define_T2(t,c,k)

%parameters
%%rate constants
kf1=k(1); kr1=k(2); k2=k(3); kf3=k(4); kr3=k(5); k4=k(6); kfold=k(7); kf5=k(8); kr5=k(9); kdeg=k(10);
kf6=k(11); kr6=k(12); k7=k(13); k8=k(14); k9=k(15); kfold2=k(16);

%%concentrations

D=c(1); DRp=c(2); Rp=c(3); mRNA=c(4); R=c(5); mRNAR=c(6); Punf=c(7); P=c(8); Xm=c(9);mRNAXm=c(10);
D2=c(11); D2P=c(12); m2=c(13); m2Xm=c(14); m2R=c(15); P2unf=c(16); P2=c(17);

%%ordinary differential equations

%Transkription
Ddot=-kf1*D*Rp+kr1*DRp+k2*DRp;
Rpdot=Ddot;
DRpdot=-Ddot;

% %Translation und Proteinfaltung
% Rdot=-kf3*mRNA*R+kr3*mRNAR+k4*mRNAR;
% mRNARdot=-Rdot;
% Punfdot=k4*mRNAR;
% %Pdot=kfold*Punf;

%leakyTranslation=Translation mit neuer Rate für kf3

%Translation und Proteinfaltung für erste und zweite mRNA
Rdot=-kf3*mRNA*R+kr3*mRNAR+k4*mRNAR-k8*m2*R+kr3*m2R+k9*m2R;
mRNARdot=kf3*mRNA*R-kr3*mRNAR-k4*mRNAR;
m2Rdot=-Rdot-mRNARdot;
Punfdot=k4*mRNAR-kfold*Punf;
P2unfdot=k9*m2R-kfold2*P2unf;
P2dot=kfold2*P2unf;

% RNA-Degradation
% Xmdot=-kf5*mRNA*Xm+kr5*mRNAXm+kdeg*mRNAXm;
% mRNAXmdot=-Xmdot;
% mRNAdot=k2*DRp+Rdot+Xmdot-kdeg*mRNAXm;

%Degradation für erste und zweiten mRNA (mRNA und m2)
Xmdot=-kf5*mRNA*Xm+kr5*mRNAXm+kdeg*mRNAXm-kf5*m2*Xm+kr5*m2Xm+kdeg*m2Xm;
mRNAXmdot=kf5*mRNA*Xm-kr5*mRNAXm-kdeg*mRNAXm;
m2Xmdot=-Xmdot-mRNAXmdot;
mRNAdot=k2*DRp-mRNARdot-mRNAXmdot-kdeg*mRNAXm;

%Transkription einer zweiten mRNA durch T7-Polymerase=P 
D2dot=-kf6*D2*P+kr6*D2P+k2*D2P;
Pdot=kfold*Punf+D2dot;
D2Pdot=-D2dot;
m2dot=k7*D2P-m2Rdot-m2Xmdot-kdeg*m2Xm;




sol=[Ddot;DRpdot;Rpdot;mRNAdot;Rdot;mRNARdot;Punfdot;Pdot;Xmdot;mRNAXmdot;D2dot;D2Pdot;m2dot;m2Xmdot;m2Rdot;P2unfdot;P2dot];