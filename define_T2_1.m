function sol=define_T2_1(t,c,k)

%parameters
%%rate constants
kf1=k(1); kr1=k(2); k2=k(3); kf3=k(4); kr3=k(5); k4=k(6); kfold=k(7); kf5=k(8); kr5=k(9); kdeg=k(10);
kf6=k(11); kr6=k(12); k7=k(13);

%%concentrations

D=c(1); DRp=c(2); Rp=c(3); mRNA=c(4); R=c(5); mRNAR=c(6); Punf=c(7); P=c(8); Xm=c(9);mRNAXm=c(10);
D2=c(11); D2P=c(12); m2=c(13);

%%ordinary differential equations

%Transkription
Ddot=-kf1*D*Rp+kr1*DRp+k2*DRp;
Rpdot=Ddot;
DRpdot=-Ddot;

%Translation und Proteinfaltung
Rdot=-kf3*mRNA*R+kr3*mRNAR+k4*mRNAR;
mRNARdot=-Rdot;
Punfdot=k4*mRNAR;
%Pdot=kfold*Punf;

%leakyTranslation=Translation mit neuer Rate für kf3

%RNA-Degradation
Xmdot=-kf5*mRNA*Xm+kr5*mRNAXm+kdeg*mRNAXm;
mRNAXmdot=-Xmdot;
mRNAdot=k2*DRp+Rdot+Xmdot-kdeg*mRNAXm;
%Degradation für erste und zweiten mRNA (mRNA und m2)

%Transkription einer zweiten mRNA durch T7-Polymerase=P 
D2dot=-kf6*D2*P+kr6*D2P+k2*D2P;
Pdot=kfold*Punf+D2dot;
D2Pdot=-D2dot;
m2dot=k7*D2P;




sol=[Ddot;DRpdot;Rpdot;mRNAdot;Rdot;mRNARdot;Punfdot;Pdot;Xmdot;mRNAXmdot;D2dot;D2Pdot;m2dot];