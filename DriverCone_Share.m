% DriverCone_Share, 8/18/20.
% Reflection function for conical bore discontinuity in duct: cylinder
% connected to divergent or convergent cone with same planar diameter.
% Calculate areas with Ware-Aki & layer peeling algorithms.

% Cylinder area is area0 (cm)
% Cone area is also area0, so a taper discontinuity only
% Missing (radial) length rb of cone.
close all, clear all
T=1/48; % sample period, ms
c=1e-3*34442; % good enough, cm/ms
Ds=c*T/2; % (radial) length increments are Ds = c * T/2 = 0.3588 cm
a0=0.4; % initial radius
area0=pi*a0^2;
Me=6; % number of plot points for cylindrical tube left of discontinuity
rrentry=(-Me:0)*Ds;
a0entry=repmat(a0,1,Me+1);
area0entry=repmat(area0,1,Me+1);
% Discontinuity in planar diameter geometry with value 1, but model RF and  
% area estimations below are based on no discontinuity in diameter.
swDiameterDiscontinuity=0; 
taperAngle=10*(pi/180); % taper angle with initial number in degrees
for swDivergent=[0,1]
  if ~swDiameterDiscontinuity
    if swDivergent % any even Nt is OK, and rb>0
      tAngle=taperAngle;
      rb=a0/sin(tAngle);
      Nt=11;
      stitl='Cylinder to Divergent cone';
    else % convergent bore, rb<0
      tAngle=-taperAngle;
      rb=a0/sin(tAngle);
      Nt=floor(abs(rb)/Ds)-1; % limit Nt so converging cone has area>0
      stitl='Cylinder to Convergent cone';
    end
    tauR=-2*rb/c; % time constant for RF
    rr=(0:Nt)*Ds; % radial distance along cone
    a=a0*(1+rr/rb); % radius of conical bore is sine function
    tt=T*(0:Nt);
    rf=(1/tauR)*exp(tt/tauR); % taper discontinuity make rf(1)~=0
    rf(1)=0.5*rf(1); % more accurate with 1/2 from sgn function at 0 = 0
    zz=rr*cos(tAngle);
  else % taper and diameter discontinuity, cylinder to cone
    swSmaller=1;
    Da=0.1; % relative change in radius (or diameter)
    if swSmaller
      Da=-Da;
    end
    aC=(1+Da)*a0; % radius of initial conical bore, +/- 10% change in diameter
    B=(a0/aC)^2; % area ratio of cylinder to cone
    if swDivergent
      rb=15;
      Nt=24;
      tauR=-(B+1)*rb/c;
      if aC>a0
        stitl='Cylinder to larger-diameter divergent cone';
      else
        stitl='Cylinder to smaller-diameter divergent cone';
      end
    else
      rb=-10;
      Nt=floor(abs(rb)/Ds); % limit Nt so that converging cone has positive area
      if rem(Nt,2)
        Nt=Nt-1; % make Nt even
      end
      tauR=-(B+1)*rb/c;
      if aC>a0
        stitl='Cylinder to larger-diameter convergent cone';
      else
        stitl='Cylinder to smaller-diameter convergent cone';
      end
    end
    rr=(0:(Nt-1))*Ds;
    mytan=tan(aC/rb);
    a=aC+rr*mytan; % radius of conical bore
    tt=T*(0:(Nt-1));
    % RF model from Martinez (1988) publications
    rf=2*B/((B+1)*tauR)*exp(tt/tauR); % taper discontinuity make rf(1)~=0
    rf(1)=0.5*rf(1)+(B-1)/(B+1)/T; % more accurate with 1/2 from sgn function
    % and includes delta function change of area at time 0
  end
  
  disp('Start layer-peeling'); 
  [rLPCy,areaLPCy]=LayerPeeling_Cylinder_Lossless(T,rf,area0,Nt+1);
  
  disp('Start Ware-Aki');
  [rWA,areaWA]=Ware_Aki(T,rf,area0,Nt+1);
  radiusWA=sqrt(areaWA/pi);
  radiusLPCy=sqrt(areaLPCy/pi);
  
  rrAll=[rrentry,rr]; % attach entry cylinder of length Me
  areaCone=pi*a.^2; % area of cone with initial value area0
  mytan=tan(tAngle);
  Dz=Ds*cos(tAngle);
  volumeCone=zeros(1,Nt+1);
  for jj=1:(Nt+1)
    ajj=a(jj);
    volumeCone(jj)=pi*Dz*(ajj^2+mytan*Dz*ajj+Dz^2*mytan^2/3);
  end
  meanareaCone=volumeCone/Dz;
  
  ms=6;
  hf11=figure(11);
  set(hf11,'DefaultLineMarkerSize',ms);
  if swDivergent
    subplot(2,2,1),
    h0=plot(zz,areaLPCy,'x',zz,areaWA,'o');
    hold on;
    Cylplot(0,zz,areaLPCy,0,0,ms,h0(1).Color); % ms has no effect here
    ha=plot(zz,areaCone,'k',rrentry,area0entry,'--k');
    ymax=ceil(max([areaLPCy,areaWA,areaCone,area0]));
    xtick1=0:1:4;
    xlim1=[-0.5,4];
    set(gca,'XLim',xlim1,'XTick',xtick1,'YLim',[0,ymax]);
    if swDiameterDiscontinuity
      plot([0,0],[area0,areaCone(1)],'k-');
    end
    text(xlim1(2),0,'A ',...
      'HorizontalAlignment','right','VerticalAlignment','bottom');
    hl=legend('Layer peeling','Ware-Aki','Cone','Cylinder',...
      'Location','NorthWest');
    set(hl,'AutoUpdate','off','box','off');
    xlabel('Axial distance (cm)');
    ylabel('Area (cm^2)');
    title({stitl,'Area'});% (geometrical and estimated)'});
    
    subplot(2,2,3),
    DareaLPCy=100*(areaLPCy-areaCone)./areaCone
    DareaWA=100*(areaWA-areaCone)./areaCone;
    DareaGeom=100*(areaLPCy-meanareaCone)./meanareaCone
    plot(zz,DareaLPCy,'x',zz,DareaWA,'o',zz,DareaGeom,'d');
    yL=[-1,20];
    set(gca,'XLim',xlim1,'XTick',xtick1,'YLim',yL,'YTick',0:5:20);
    text(xlim1(2),yL(2),'B ',...
      'HorizontalAlignment','right','VerticalAlignment','top');
    hl=legend('Layer peeling','Ware-Aki','Layer peeling re: mean area',...
      'Location','NorthEast');
    set(hl,'AutoUpdate','off','box','off');
    xlabel('Axial distance (cm)');
    ylabel('Rel. Error (%)');
    title({stitl,'Relative Error in Area'});
  else
    subplot(2,2,2),
    h0=plot(zz,areaLPCy,'x',zz,areaWA,'o');
    hold on;
    Cylplot(0,zz,areaLPCy,0,0,ms,h0(1).Color)
    ha=plot(zz,areaCone,'k',rrentry,area0entry,'--k');
    ymax=ceil(max([areaLPCy,areaWA,areaCone,area0]));
    xlim2=[-0.5,2];
    xtick2=0:0.5:2;
    set(gca,'XLim',xlim2,'XTick',xtick2,...
      'YLim',[-0.1,ymax],'YTick',0:0.25:1);
    text(xlim2(2),ymax,'C ',...
      'HorizontalAlignment','right','VerticalAlignment','top');
    if swDiameterDiscontinuity
      plot([0,0],[area0,areaCone(1)],'k-');
    end
    xlabel('Axial distance (cm)');
    ylabel('Area (cm^2)');
    title({stitl,'Area'});% (geometrical and estimated)'});
    
    subplot(2,2,4),
    DareaLPCy=100*(areaLPCy-areaCone)./areaCone
    DareaWA=100*(areaWA-areaCone)./areaCone;
    DareaGeom=100*(areaLPCy-meanareaCone)./meanareaCone
    plot(zz,DareaLPCy,'x',zz,DareaWA,'o',zz,DareaGeom,'d');
    yL=[-80,5];
    set(gca,'XLim',xlim2,'XTick',xtick2,'YLim',yL);%,'YTick',-20:-20:-80);
    text(xlim2(2),yL(2),'D ',...
      'HorizontalAlignment','right','VerticalAlignment','top');
    xlabel('Axial distance (cm)');
    ylabel('Rel. Error (%)');
    title({stitl,'Relative Error in Area'});
  end
  
end

