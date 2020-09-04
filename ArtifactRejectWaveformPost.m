function iRng = ArtifactRejectWaveformPost(trials)
% ArtifactRejectWaveformPost. 3/28/2019. 
% Used in calibration stage to reject tube outliers based on peak test.
%Constuct median of sum over time 
% of the absolute deviations of the median, and then perform a MAD test
% over that across the individual trials to identify bad trials.  iRng is
% the good trials.
MADoutlier=1.5; % 3 is Tukey "near outlier" value
% Calculate outliers from median absolute deviation (MAD) test.
MADwaveform=max(trials,[],2);
MAD=median(MADwaveform);
fenceMAD=MADoutlier*median(abs(MADwaveform-MAD));
iRng=abs(MADwaveform-MAD)<=fenceMAD; % iRng of good trials
end
