function totalError=ode23s_solver_T2(param,c_0,data,errW,cutFFT)

global bestparam besterr
k=param;

deadTime=1; %dataPoints
deltaT=3*60; % in sec
endTime=22; %dataPoints 250 bis zum ende
tSpanSim= [deadTime*deltaT:deltaT:endTime*deltaT];

expData=data(deadTime:endTime)';

options=odeset('AbsTol',1e-9); % set tolerances

[simTime,simData]=ode23s(@(t,c)define_T2(t,c,k),tSpanSim,c_0,options);

simData1=simData(:,1);
simData2=simData(:,2);
simData3=simData(:,3);
simData4=simData(:,4);
simData5=simData(:,5);
simData6=simData(:,6);
simData7=simData(:,7);
simData8=simData(:,8);
simData9=simData(:,9);
simData10=simData(:,10);
simData11=simData(:,11);
simData12=simData(:,12);
simData13=simData(:,13);
simData14=simData(:,14);
simData15=simData(:,15);
simData16=simData(:,16);
simData17=simData(:,17);


plot(simTime,simData1,simTime,simData2,simTime, simData3, simTime, simData4,simTime, simData5,...
    simTime, simData6,simTime, simData7,simTime, simData8,simTime, simData9,simTime, simData10,...
    simTime, simData11,simTime, simData12,simTime, simData13,simTime, simData14,simTime, simData15,...
    simTime, simData16,simTime, simData17,simTime,expData)
legend('D','DRp','Rp','mRNA','R','mRNAR','Punf','P','Xm','mRNAXm','D2','D2P','m2','m2Xm','m2R','P2unf','P2','expData')



%% 
von=1; bis=length(expData); l1=2000; 
FFTregion=zeros(l1,1);
FFTregion(von:bis,1)=1;

% readout of simulated data
expData=expData/max(expData(find(FFTregion))); 

simDataP2=simData17/max(simData17(find(FFTregion)));

expDataInFFT=expData(find(FFTregion))/max(expData(find(FFTregion)));
simDataInFFT=simDataP2(find(FFTregion))/max(simDataP2(find(FFTregion)));


simTimeInFFT=simTime(find(FFTregion)); %usually simTime should be the same as expTime
expTimeInFFT=simTimeInFFT;

%% 
lsqE=lsqError(expData,simData17);
[fftE,dExpdT,dSimdT,tExpFFT,tSimFFT] = fftError(expDataInFFT,simDataInFFT,expTimeInFFT,simTimeInFFT,cutFFT);

%weight factors for errors
lsqW=errW(1); fftW=errW(2);

totalError=lsqE*lsqW+fftE*fftW;

%% update error and parameters if total error is reduced
if totalError<besterr
    bestparam=param;
    besterr=totalError
end

