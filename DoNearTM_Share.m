% DoNearTM_Share, 5/28/20 based on DoNearTM, 2/5/20.  
% Driver for area-distance function estimation and sound field near TM 
%   based on measured RF at probe tip.

close all
clear all
% Load data files
% Calibration params: rel level 10 dB,nMat=72,nPmultiple=4,swChi1=1,swEvanes=0,CGLS
fnDataRoot='Cal-20190502T124516_Out_20190515T184842'; % RF calibration file
load(fnDataRoot,'p0NormInc','ipMax1','pMax1','nMat2p','tmsP0');
% Ear data file
fnEarRoot='DhkBlue-p10_R_A_more'; % choice for results in 2020 JASA paper
load(fnEarRoot,'uprm','xpr','RFear','Naf');
% Initial processing
Base=BaseRFShare(0);
Base.TubeArea=pi*0.4^2; % area of tubes used in calibration, not test system
Base.temperature=uprm.Info.Temperature;
Base.altitude=uprm.Info.Altitude; % reset Base.altitude=0 for sea-level output
Base.nnE=xpr.Parameters.BufferSize;
sTitl=strrep(fnEarRoot,'_','\_');
% Calc loss-less area functions by Ware-Aki and cylindrical layer peeling,
%   and use to calculate Mnear
RF_TubeArea=RFear;
MnearMax=10; % max number of samples for initial area function methods
rf=RF_TubeArea(1:MnearMax);
[rLPCy,areaLPCy]=LayerPeeling_Cylinder_Lossless(Base.T,rf,Base.TubeArea,...
  MnearMax);
[rWA,areaWA]=Ware_Aki(Base.T,rf,Base.TubeArea,MnearMax);
dArea=(areaWA-areaLPCy);
mendTM=find(areaLPCy<0.1,1,'first'); %1st one at small end
if mendTM>1
  Mnear=mendTM-2; 
end
if isempty(Mnear) || Mnear<2
  disp('Error condition');
  return
end
hRF=Base.T*RFear; % dimensionless RF at probe tip
ifig=1;
hf0=figure(ifig);
ifig=ifig+1;
co=colororder;
for jj=1:2
  swLossy=jj-1; % 0 for jj=1 and 1 for jj==2
  if swLossy
    ls='-';
  else
    ls='--';
  end
  % Forward transfer function to near TM
  [grP,grM,glP,glM,areaD,rD,nOffset]=CalcNearTM(Base,swLossy,Naf,hRF,...
    MnearMax);
  rhoc=1; % units scaled by 1/rhoc
  Yc0=1;
  Yc=areaD/Base.TubeArea/rhoc; % Yc is dimensionless here
  Mm1=Mnear-1;
  ixMm1=transpose(1:Mm1);
  YcN=Yc(ixMm1);
  WP=[Yc0;YcN.*sum(grP(ixMm1,:).^2,2)];
  WM=[-Yc0*sum(hRF.^2);-YcN.*sum(grM(ixMm1,:).^2,2)];
  W=WP+WM;
  xx=[1;ixMm1+1]; % right segment m-1 is at node m
  xlim=[1,Mm1+1];
  hp=plot(xx,WP,'o',xx,-WM,'x',xx,W,'d');
  for jj=1:length(hp)
    hp(jj).Color=co(jj,:);
    hp(jj).LineStyle=ls;
  end
  xlabel('Node number m');
  ylabel('Normalized Energy');
  if swLossy
    legend('Forward','-1 \times Reverse','Absorbed',...
      'Location','SouthEast','AutoUpdate','off','Box','off');
  end
  set(gca,'XLim',[0.8,Mnear+0.2],'XTick',xx,...
    'YLim',[0,1.05],'YTick',0:0.2:1);
  drawnow;
  h=findall(hf0, '-property', 'LineWidth');
  lw=1.2;
  set(h,{'LineWidth'}, num2cell(lw))
  h=findall(hf0, '-property', 'fontsize');
  hFontSize = cell2mat(get(h,'FontSize'));
  newFontSize =hFontSize *1.2;
  set(h,{'FontSize'}, num2cell(newFontSize));
  drawnow;
  hold on;

  ifig=PlotNearTM_Share(ifig,Base,swLossy,Naf,hRF,grP,grM,glP,Mnear,...
    nOffset,tmsP0,nMat2p,pMax1,p0NormInc);
end

