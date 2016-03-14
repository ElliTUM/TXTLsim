function plot3Dprofil

IPTG=[0:50000:600000];
aTc=[0:50000:600000];

for i=1:length(IPTG)
    for j=1:length(aTc)
    I=IPTG(i);
    I2=aTc(j);
   main_Komp_2(I,I2);
   mCer(i,j)= ans(14);
    end
  
end

%% create colormap
n=50;
for i=1:n
   
     ColorI2{i}=rgb('black') - (i-1)*(rgb('black')-rgb('green'))/n;
end
for i1=1:n
    for i2=1:3
     
        map2(i1,i2)=ColorI2{1,i1}(1,i2);
    end
end
%%
LWidth=3; % Linewidth
FontForPlot='Arial'; % Font type
FSize=24; % Font size
MSize=3; % Marker size
%%
figure
colormap(map2)
surf(IPTG,aTc,mCer,'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','gouraud','Linewidth',LWidth)
xlim([400 1600])
axis tight
title('mCerulean(IPTG,aTc)')
xlabel('IPTG (nM)')
ylabel('aTC (nM)')
zlabel('mCerulean')
view(-10,30)
%camlight left
%camlight right
brighten(0.65)
set(gca,'FontName',FontForPlot,'FontSize',FSize)