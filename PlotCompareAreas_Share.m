function ifig=PlotCompareAreas_Share(ifig,TubeArea,z,zTM,zendTM,...
  areaLPCy,dArea,dAreaLossy,sTitl)
% Plot outputs from cylindrical layer-peeling and Ware-Aki methods 

xlim=[-1,z(end)+1];
xtick=0:3:(3*floor(1/3*xlim(2)));
ms=4;
hf=figure(ifig);
ifig=ifig+1;
set(hf,'Position',[200,200,420,680]);
subplot(3,1,1),  
Cylplot(1e2*TubeArea,z,1e2*areaLPCy,xlim,xtick,ms);
set(gca,'XLim',xlim,'XTick',xtick);
hold on;
if length(sTitl)>=4
  if strcmp(sTitl(1:4),'Tube')
    plot(xlim,[EarArea,EarArea],'k:');
    title(['Geometrical area of tube: ',num2str(EarArea,'%6.3f'),' cm^2']);
  end
end
xlabel('Distance (mm)');
ylabel('Area  (mm^2)');
ylim=get(gca,'YLim');
fs=9;
text(0,0.5*ylim(2),{' Probe-tip of canal',' (m=1)'},'FontSize',fs,...
  'HorizontalAlignment','left','VerticalAlignment','top');
text(zTM,0.97*ylim(2),'(m=6) ','FontSize',fs,...
  'HorizontalAlignment','right','VerticalAlignment','top');
text(zTM,0.97*ylim(2),{' Near-TM',' region'},'FontSize',fs,...
  'HorizontalAlignment','left','VerticalAlignment','top');
plot([zTM,zTM],ylim,'k:');
plot([zendTM,zendTM],ylim,'k:');
zEndFreq=[25.4,25.4];
plot(zEndFreq,[26,46],'k--');
title('Plane-wave layer-peeling (loss-less)');
text(zEndFreq(1),36,{' Length',' estimate'},'FontSize',fs,...
  'HorizontalAlignment','left','VerticalAlignment','middle');
dxt=xlim(2)-0.2;
text(dxt,ylim(2),'A',...
  'HorizontalAlignment','right','VerticalAlignment','top');

text(zendTM,0.41*ylim(2),{'  End-TM','  region (m=8)'},'FontSize',fs,...
  'HorizontalAlignment','left','VerticalAlignment','top');
subplot(3,1,2),
s1='B';
Cylplot(0,z,1e2*dArea,xlim,xtick,ms);
hold on;
plot(xlim,[0,0],'k:');
ylimd=[-0.7,0.5];
plot([zTM,zTM],ylimd,'k:');
plot([zendTM,zendTM],ylimd,'k:');
ytickd=-0.6:0.2:0.4;
set(gca,'XLim',xlim,'XTick',xtick,'YLim',ylimd,'YTick',ytickd);
title('Ware-Aki minus Layer-peeling');
xlabel('Distance (mm)');
ylabel('\Delta Area  (mm^2)');
text(dxt,ylimd(2),s1,...
  'HorizontalAlignment','right','VerticalAlignment','top');

subplot(3,1,3),
s1='C';
Cylplot(0,z,1e2*dAreaLossy,xlim,xtick,ms);
hold on;
plot(xlim,[0,0],'k:');
plot([zTM,zTM],ylimd,'k:');
plot([zendTM,zendTM],ylimd,'k:');
set(gca,'XLim',xlim,'XTick',xtick,'YLim',ylimd,'YTick',ytickd);
title('Lossy minus Loss-less');
xlabel('Distance (mm)');
ylabel('\Delta Area  (mm^2)');
text(dxt,ylimd(2),s1,...
  'HorizontalAlignment','right','VerticalAlignment','top');

  function Cylplot(TubeArea,z,area,xlim,xtick,ms)
    z0=[xlim(1),xtick(1)];
    plot(z0,repmat(TubeArea,1,2),':k');
    hold on;
    N=length(z);
    plot(repmat(z(1),1,2),[TubeArea,area(1)],':k');
    for kk=1:(N-1)
      dz=[z(kk),z(kk+1)];
      plot(dz,repmat(area(kk),1,2),'k-',...
        repmat(z(kk+1),1,2),[area(kk),area(kk+1)],'k-');
      plot(z(kk),area(kk),'ko','MarkerSize',ms);
    end
    plot(z(N),area(N),'ko','MarkerSize',ms);
  end
end