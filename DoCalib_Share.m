% DoCalib_Share, 8/18/20.  Script to call LoadProcessTubeData_RF_Share to 
% plot various model reflection functions for a long and short cylindrical
% tube.
close all
set(0,'DefaultLineLineWidth',0.5);
fnDataRoot='Cal-20190502T124516';%Rel level 10 dB, H stm, for Keefe (2020b)
load(fnDataRoot,'TubeData','TubeList','uprm','xpr');
% Structures written out by software for 2020b paper
% uprm, user parameters
% xpr, experiment list parameters
% TubeList, tube list parameters
% TubeData, chirp data acquired in calibration tubes
Base=BaseRFShare(0); % Base object
Base.temperature=uprm.Info.Temperature;
Base.altitude=uprm.Info.Altitude;% reset Base.altitude=0 for sea-level output
LoadProcessTubeData_RF_Share(Base,TubeData,TubeList,xpr);
