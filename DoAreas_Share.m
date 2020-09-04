% DoAreas_Share 8/18/20.  Driver for area-distance function estimation 
%   and estimate sound field near the TM.

close all
clear all
% Load data files
% Calibration params: rel level 10 dB,nMat=72,nPmultiple=4,swChi1=1,swEvanes=0,CGLS
fnDataRoot='Cal-20190502T124516_Out_20190515T184842'; % RF calibration file
load(fnDataRoot,'nMat2p','tmsP0');
% Ear data file
fnEarRoot='DhkBlue-p10_R_A_more'; % choice for results in 2020 JASA paper
load(fnEarRoot,'uprm','xpr','RFear','Naf');
% Initial processing
Base=BaseRFShare(0); % for viscothermal and general run parameter values
Base.TubeArea=pi*0.4^2; % area of tubes used in calibration, not test system
Base.temperature=uprm.Info.Temperature;
Base.altitude=uprm.Info.Altitude;  % reset Base.altitude=0 for sea-level output
sTitl=strrep(fnEarRoot,'_','\_');
% Calc loss-less area functions by Ware-Aki and cylindrical layer peeling
% Used to calculate Mnear
RF_TubeArea=RFear;
MnearMax=10; % max number of samples for initial area function methods
rf=RF_TubeArea(1:MnearMax);
[rLPCy,areaLPCy]=LayerPeeling_Cylinder_Lossless(Base.T,rf,...
  Base.TubeArea,MnearMax); % from Sharp dissertation
[rWA,areaWA]=Ware_Aki(Base.T,rf,Base.TubeArea,MnearMax);
dArea=(areaWA-areaLPCy);
mendTM=find(areaLPCy<0.1,1,'first'); %1st one at small end
if mendTM>1
  Mnear=mendTM-2; % 2 Ds=7.2 mm about the spatial extent of the TM
end
if isempty(Mnear) || Mnear<2
  disp('Error condition');
  return
end
% Forward transfer function to near TM, only to estimate lossy areaD
[~,~,~,~,areaD]=CalcNearTM(Base,1,Naf,Base.T*RFear,MnearMax);

areaNoLoss=transpose(areaLPCy);
WareAki=transpose(areaWA);
dAreaLossy=areaD-areaNoLoss;
Mnear
Ds=10*Base.T*Base.airVec.c/2; % Ds in mm
z=(0:(MnearMax-1))*Ds;
zTM=z(Mnear)+0.15*Ds;
zendTM=z(mendTM)-0.15*Ds;
set(0,'DefaultLegendBox','off','DefaultLegendAutoUpdate','off');
ifig=PlotCompareAreas_Share(1,Base.TubeArea,z,zTM,zendTM,areaLPCy,...
  dArea,dAreaLossy,sTitl);
