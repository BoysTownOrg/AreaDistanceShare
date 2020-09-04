function ifig=PlotNearTM_Share(ifig,Base,swLossy,Naf,hRF,grP,grM,glP,Mnear,...
  nOffset,tmsP0,nMat2p,pMax1,p0NormInc)
% CalcNearTM_Share, 5/28/2020 based on CalcNearTM, 2/26/2020.
% Calculate area-distance function from measured RF
% at probe tip (hRF) over duration of Naf samples for discontinuities from
% 1:MnearMax with or without a model of viscothermal losses (swLossy).
% data acquired by AcquireData_RFv0

T=Base.T;
Tms=Base.Tms;
MTMf=Mnear;
MTMr=Mnear-1;
ImpulseF=zeros(1,Naf);
ImpulseF(1)=1;
gTotTM=grP(MTMr,:)+grM(MTMr,:);
swLevel=1; %   used within CalcSpec
domegakHz=2*pi*Base.fskHz/Base.mMat;  % used within CalcSpec
[LFp0,GDFp0]=CalcSpec(ImpulseF,1); % forward incident impulse
[LFp,GDFp]=CalcSpec(glP(1,:),1); % glP  for TF out from probe tip
[LFtm,GDFtm]=CalcSpec(glP(MTMf,:),1); % glP  for TF into near-TM
[LRp,GDRp]=CalcSpec(hRF,1);
[LRtm,GDRtm]=CalcSpec(grM(MTMr,:),1);
[LTp,GDTp]=CalcSpec(ImpulseF+hRF,1);
[LTtm,GDTtm]=CalcSpec(gTotTM,1);
[LFtmM1,~]=CalcSpec(grP(MTMr,:),1); %incoming forward wave to near-TM
MTMdelay=nOffset(MTMr)*Tms;
GDoffset=repmat(MTMdelay,1,length(GDFtm));
GDFtm=GDFtm+GDoffset;
GDRtm=GDRtm+GDoffset;
GDTtm=GDTtm+GDoffset;
mydir=[fileparts(which(mfilename)),filesep];
fn=[mydir,'LosslessResults'];
if swLossy==1
  s=load(fn,'GDFtm','GDRtm','LFtm','LRtm','GDFp','GDRp','LFp','LRp');
else
  save(fn,'GDFtm','GDRtm','LFtm','LRtm','GDFp','GDRp','LFp','LRp');
  return % no plots made for loss-less
end

x=(0:(Naf-1))*Tms;
xTM=nOffset(Mnear)*Tms+x;
xlim=[-3,60]*Tms;
xtick=0:0.2:1.2;
ylim=[-0.3,1.3];
ytick=-0.2:0.2:1.2;
xOffset=[0,nOffset(Mnear),nOffset(Mnear)]*Tms;
gzOffset=[0,0,1];
co=get(groot,'defaultAxesColorOrder');
cb=co(1,:);
cr=co(2,:);

hf1=figure(ifig);
ifig=ifig+1;
set(hf1,'Position',[100,100,360,760]);
subplot(3,1,3)
hp=plot(x,ImpulseF+hRF,xTM,gTotTM,xOffset,gzOffset*gTotTM(1));
hp(1).Color=cr;
hp(2).Color=cb;
hp(3).Color=cb;
xlabel('Time  (ms)');
ylabel('Amplitude');
title('Total pressure');
set(gca,'XLim',xlim,'XTick',xtick,'YLim',ylim,'YTick',ytick);
legend('Probe tip','Near TM','Location','North','AutoUpdate','off','Box','off');
dx8=0.02;
text(xlim(2)-dx8,ylim(2),'C',...
  'HorizontalAlignment','right','VerticalAlignment','top');

subplot(3,1,1)
hp=plot(x-0.001,ImpulseF,'--',x,glP(1,:),'-',xTM,glP(MTMf,:),'-',...
  xOffset,gzOffset*glP(MTMf,1),'-');
hp(1).Color='k';
hp(2).Color=cr;
hp(3).Color=cb;
hp(4).Color=cb;
set(gca,'XLim',xlim,'XTick',xtick,'YLim',ylim,'YTick',ytick);
title('Forward pressure');
xlabel('Time  (ms)');
ylabel('Amplitude');
legend('Forward impulse probe tip','Probe tip','Near TM',...
  'Location','North','AutoUpdate','off','Box','off');
text(xlim(2)-dx8,ylim(2),'A',...
  'HorizontalAlignment','right','VerticalAlignment','top');

subplot(3,1,2)
hp=plot(x,hRF,xTM,grM(MTMr,:),xOffset,gzOffset*grM(MTMr,1));
hp(1).Color=cr;
hp(2).Color=cb;
hp(3).Color=cb;
set(gca,'XLim',xlim,'XTick',xtick,'YLim',ylim,'YTick',ytick);
title('Reverse pressure');
xlabel('Time  (ms)');
ylabel('Amplitude');
legend('Probe tip','Near TM','Location','North','AutoUpdate','off','Box','off');
text(xlim(2)-dx8,ylim(2),'B',...
  'HorizontalAlignment','right','VerticalAlignment','top');

logfr=log2(1e-3*Base.fAnalysis); % for plotting re:0 at 1 kHz
XLimf=[-2, 4];
XTickf=-2:4;
XTickLbl={'0.25','0.5','1','2','4','8','16'};
if swLevel
  sy='Level (dB)';
else
  sy='';
end

hf2=figure(ifig);
ifig=ifig+1;
set(hf2,'Position',[100,100,760,840]);

ixz=find(LRp>0,1);
subplot(3,2,5)
hp=plot(logfr,LTp,logfr,LTtm);
hp(1).Color=cr;
hp(2).Color=cb;
if swLevel
  ylimLf=[-12,7.5];
  ytickLf=-12:3:6;
  dx9=XLimf(1)+0.2;
else
  ylimLf=[0,3];
  ytickLf=0:0.5:3;
  dx9=XLimf(1)+0.1; % for magnitude
end
set(gca,'XLim',XLimf,'XTick',XTickf,'XTickLabel',XTickLbl,...
  'YLim',ylimLf,'YTick',ytickLf);
title('Total pressure level');
xlabel('Frequency (kHz)');
ylabel(sy);
text(dx9,ylimLf(2),'E',...
  'HorizontalAlignment','left','VerticalAlignment','top');
legend('Probe tip','Near TM','Location','SouthWest','AutoUpdate','off','Box','off');
hold on;
plot(repmat(logfr(ixz),1,2),ylimLf,'r:');
plot(XLimf,repmat(6,1,2),'k:');

subplot(3,2,6), % add delay of xTM in time
hp=plot(logfr,GDTp,'-',logfr,GDTtm,'-',XLimf,repmat(MTMdelay,1,2),'k:');
hp(1).Color=cr;
hp(2).Color=cb;
ylimgd=[-0.5,0.5];
ytickgd=-0.4:0.2:0.4;
set(gca,'XLim',XLimf,'XTick',XTickf,'XTickLabel',XTickLbl,...
  'YLim',ylimgd,'YTick',ytickgd);
title('Total pressure group delay');
xlabel('Frequency (kHz)');
ylabel('Group delay (ms)');
text(dx9,ylimgd(2),'F',...
  'HorizontalAlignment','left','VerticalAlignment','top');
hl=legend('Probe tip','Near TM','One-way delay to Near TM',...
  'Location','SouthWest','AutoUpdate','off','Box','off');
hold on;
plot(repmat(logfr(ixz),1,2),ylimLf,'r:');

subplot(3,2,3)
hp=plot(logfr,LRp,logfr,LRtm);
hp(1).Color=cr;
hp(2).Color=cb;
set(gca,'XLim',XLimf,'XTick',XTickf,'XTickLabel',XTickLbl,...
  'YLim',ylimLf,'YTick',ytickLf);
title('Reverse pressure level');
xlabel('Frequency (kHz)');
ylabel(sy);
text(dx9,ylimLf(2),'C',...
  'HorizontalAlignment','left','VerticalAlignment','top');
hold on;
plot(logfr,LFtmM1,'k-.');
legend('Probe tip (measured RF)','Near TM',...
  'Forward pressure into Near TM',...
  'Location','SouthWest','AutoUpdate','off','Box','off');
plot(XLimf,zeros(1,2),'k:');
plot(repmat(logfr(ixz),1,2),ylimLf,'r:');
text(logfr(ixz),6,['Probe-tip level>0 dB at ',...
  num2str(2^logfr(ixz),'%.2f'),' kHz'],'HorizontalAlignment','right');

subplot(3,2,4),
hp=plot(logfr,GDRp,'-',logfr,GDRtm,'-',...
  logfr,GDRp-repmat(MTMdelay,1,length(logfr)),'k-.');
hp(1).Color=cr;
hp(2).Color=cb;
set(gca,'XLim',XLimf,'XTick',XTickf,'XTickLabel',XTickLbl,...
  'YLim',ylimgd,'YTick',ytickgd);
title('Reverse pressure group delay');
xlabel('Frequency (kHz)');
ylabel('Group delay (ms)');
text(dx9,ylimgd(2),'D',...
  'HorizontalAlignment','left','VerticalAlignment','top');
legend('Probe tip (measured RF)','Near TM',...
  'Probe tip minus one-way delay to Near TM','Location','SouthWest','AutoUpdate','off','Box','off');
hold on;
plot(repmat(logfr(ixz),1,2),ylimLf,'r:');

subplot(3,2,1)
hp=plot(logfr,LFp0,'-.',logfr,LFp,'-',logfr,LFtm,'-');
hp(1).Color='k';
hp(2).Color=cr;
hp(3).Color=cb;
set(gca,'XLim',XLimf,'XTick',XTickf,'XTickLabel',XTickLbl,...
  'YLim',ylimLf,'YTick',ytickLf);
title('Forward pressure level');
xlabel('Frequency (kHz)');
ylabel(sy);
text(dx9,ylimLf(2),'A','HorizontalAlignment','left',...
  'VerticalAlignment','top');
legend('Incident','Probe tip','Near TM',...
  'Location','SouthWest','AutoUpdate','off','Box','off');
hold on;
plot(repmat(logfr(ixz),1,2),ylimLf,'r:');
plot([logfr(ixz)-1/3,logfr(ixz)+1/3],repmat(LFtm(ixz),1,2),'k:',...
  logfr(ixz),LFtm(ixz),'ko');
text(logfr(ixz)-1/4,LFtm(ixz),[num2str(LFtm(ixz),'%.1f'),' dB'],...
  'HorizontalAlignment','right','VerticalAlignment','top');
[~,ixmax]=max(LFtm);
plot([logfr(ixmax)-1/3,logfr(ixmax)+1/3],repmat(LFtm(ixmax),1,2),'k:',...
  logfr(ixmax),LFtm(ixmax),'ko');
text(logfr(ixmax),LFtm(ixmax),[num2str(LFtm(ixmax),'%.1f'),' dB'],...
  'HorizontalAlignment','center','VerticalAlignment','bottom');
hold on;

subplot(3,2,2),
hp=plot(logfr,GDFp0,'k-.',logfr,GDFp,'k-',logfr,GDFtm,'k-',...
  XLimf,repmat(MTMdelay,1,2),'k:');
hp(1).Color='k';
hp(2).Color=cr;
hp(3).Color=cb;
set(gca,'XLim',XLimf,'XTick',XTickf,'XTickLabel',XTickLbl,...
  'YLim',ylimgd,'YTick',ytickgd);
title('Forward pressure group delay');
xlabel('Frequency (kHz)');
ylabel('Group delay (ms)');
text(dx9,ylimgd(2),'B',...
  'HorizontalAlignment','left','VerticalAlignment','top');
hl=legend('Incident','Probe tip','Near TM',...
  'One-way delay to Near TM',...
  'Location','SouthWest','AutoUpdate','off','Box','off');
hold on;
plot(repmat(logfr(ixz),1,2),ylimLf,'r:');

rngp0=1:min(length(tmsP0),nMat2p); % needed for area function
p0mPa=pMax1*p0NormInc(rngp0); % incident click tube 1 pressure in mPa
pEarmPaF=1e-3*conv(p0mPa,glP(MTMf,:));
pEarmPaR=1e-3*conv(p0mPa,grM(MTMr,:));
Nc=length(pEarmPaF);
xc=0:(Nc-1);
xc=xc*Tms;
xlimc=[0,Nc*Tms];

% Increase linewidth and font size
drawnow;
h=findall(hf2, '-property', 'LineWidth');
lw=1.2;
set(h,{'LineWidth'}, num2cell(lw))
drawnow;
h=findall(hf2, '-property', 'fontsize');
hFontSize = cell2mat(get(h,'FontSize'));
newFontSize =hFontSize *1.1;
set(h,{'FontSize'}, num2cell(newFontSize));
if swLossy
  hf5=figure(ifig);
  ifig=ifig+1;
  set(hf5,'Position',[105,105,760,560]);
  subplot(2,2,1),
  hp=plot(logfr,LFp-s.LFp,logfr,LFtm-s.LFtm);
  hp(1).Color=cr;
  hp(2).Color=cb;
  ylimfLD=[-0.5,0.5];
  ytickfLD=-0.4:0.2:0.4;
  set(gca,'XLim',XLimf,'XTick',XTickf,'XTickLabel',XTickLbl,...
    'YLim',ylimfLD,'YTick',ytickfLD);
  title(['\Delta','Forward pressure level']);
  xlabel('Frequency (kHz)');
  ylabel(['\Delta',sy]);
  text(dx9,0.5,'A',...
    'HorizontalAlignment','left','VerticalAlignment','top');
  legend('Probe tip','Near TM',...
    'Location','SouthWest','AutoUpdate','off','Box','off');
  hold on;
  plot(repmat(logfr(ixz),1,2),ylimLf,'r:');
  
  ylimgd=[-0.0025,0.0025];
  ytickgd=-0.002:0.001:0.002;
  YTickLbl={'-0.002','-0.001','0','0.001','0.002'};
  subplot(2,2,2),
  hp=plot(logfr,GDFp-s.GDFp,'-',logfr,GDFtm-s.GDFtm,'-');
  hp(1).Color=cr;
  hp(2).Color=cb;
  set(gca,'XLim',XLimf,'XTick',XTickf,'XTickLabel',XTickLbl,...
    'YLim',ylimgd,'YTick',ytickgd,'YTickLabel',YTickLbl);
  title(['\Delta','Forward pressure group delay']);
  xlabel('Frequency (kHz)');
  ylabel(['\Delta','Group delay (ms)']);
  legend('Probe tip','Near TM',...
    'Location','SouthWest','AutoUpdate','off','Box','off');
  hold on;
  plot(repmat(logfr(ixz),1,2),ylimLf,'r:');
  text(dx9,ylimgd(2),'B',...
    'HorizontalAlignment','left','VerticalAlignment','top');
  
  subplot(2,2,3),
  hp=plot(logfr,LRp-s.LRp,logfr,LRtm-s.LRtm);
  hp(1).Color=cr;
  hp(2).Color=cb;
  set(gca,'XLim',XLimf,'XTick',XTickf,'XTickLabel',XTickLbl,...
    'YLim',ylimfLD,'YTick',ytickfLD);
  title(['\Delta','Reverse pressure level']);
  xlabel('Frequency (kHz)');
  ylabel(['\Delta',sy]);
  text(dx9,0.5,'C',...
    'HorizontalAlignment','left','VerticalAlignment','top');
  legend('Probe tip (measured RF)','Near TM',...
    'Location','SouthWest','AutoUpdate','off','Box','off');
  hold on;
  plot(repmat(logfr(ixz),1,2),ylimLf,'r:');
  
  subplot(2,2,4),
  hp=plot(logfr,GDRp-s.GDRp,'-',logfr,GDRtm-s.GDRtm,'-');
  hp(1).Color=cr;
  hp(2).Color=cb;
  set(gca,'XLim',XLimf,'XTick',XTickf,'XTickLabel',XTickLbl,...
    'YLim',ylimgd,'YTick',ytickgd,'YTickLabel',YTickLbl);
  title(['\Delta','Reverse pressure group delay']);
  xlabel('Frequency (kHz)');
  ylabel(['\Delta','Group delay (ms)']);
  legend('Probe tip (measured RF)','Near TM',...
    'Location','SouthWest','AutoUpdate','off','Box','off');
  hold on;
  plot(repmat(logfr(ixz),1,2),ylimLf,'r:');
  text(dx9,ylimgd(2),'D',...
    'HorizontalAlignment','left','VerticalAlignment','top');
  
  drawnow;
  h=findall(hf5, '-property', 'LineWidth');
  lw=1.2;
  set(h,{'LineWidth'}, num2cell(lw))
  drawnow;
  h=findall(hf5, '-property', 'fontsize');
  hFontSize = cell2mat(get(h,'FontSize'));
  newFontSize =hFontSize *1.1;
  set(h,{'FontSize'}, num2cell(newFontSize));
  drawnow;
end

  function [LNorm,GD]=CalcSpec(pt,denom)
    % code assumes that N is even
    p=zeros(1,Base.mMat);
    p(1:length(pt))=pt;
    pDFT=fft(p);
    pNormF=pDFT(Base.ixfAnalysis)/denom;
    if swLevel % level in dB
      LNorm=10*log10(abs(pNormF).^2);
    else % linear squared magnitude
      LNorm=abs(pNormF).^2;
    end
    GD=-gradient(unwrap(angle(pNormF)))/domegakHz;
  end
end