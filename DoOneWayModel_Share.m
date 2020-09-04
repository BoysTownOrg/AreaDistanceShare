% DoOneWayModel, 2/19/20.  Driver to calculate and plot one-way
% viscothermal transmission in 1 m tube.

close all
clear all
ifig=1;
% Load ear data file
fnEarRoot='DhkBlue-p10_R_A_more'; % choice for results in 2020b JASA paper
load(fnEarRoot,'uprm','xpr','RFear','Naf');

Base=BaseRFShare(0);
Base.TubeArea=pi*0.4^2; % area of calibration tube
Base.temperature=uprm.Info.Temperature; 
Base.altitude=uprm.Info.Altitude; % reset Base.altitude=0 for sea-level output
Base.nnE=xpr.Parameters.BufferSize;
T=Base.T; % 1/fs;
MnearMax=279; % max number of samples for initial area function methods, 1 m length

% Forward-anechoic transfer function to near TM
Length=(MnearMax-1)*Base.airVec.c*T/2; % in cm
sL=['Length ',sprintf('%4.2f',Length),' cm'];
Naf=2*MnearMax+2;
glP=zeros(MnearMax,Naf); % g_l^+
grP=zeros(MnearMax,Naf); % g_r^+[m,n] at time 1/2 sample later than g_l^+[m,n]
%  so g_r^+[m,n] represents g_r^+[m,n+1/2]
areaCurrent=Base.TubeArea;
Nmax=Naf;
hf1=figure(1);
LossyPropagation(Base,MnearMax,sL,glP,grP,areaCurrent);

