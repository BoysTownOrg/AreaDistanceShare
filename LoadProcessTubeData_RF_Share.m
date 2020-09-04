function [pNorm,estTubeLengths,MRT,ifig]=...
  LoadProcessTubeData_RF_Share(Base,TubeData,TubeList,xpr)
% Calculate equivalent click pressure pNorm in time domain from chirp
% responses measured in four calibration tubes in Keefe (2020b) with 
% probe microphone after completing calibration to enable measurement of a 
% reflection function of the ear. The MRT holds model RF's of the four 
% calibration tubes, of which only the final tube Base.NTubes==4 is 
% relevant here. RF plots are generated for a  
Base.nnE=xpr.Parameters.BufferSize;
nTrials=xpr.Parameters.NumberOfTrials;
nBuffers=xpr.Parameters.NumberOfMeanBuffers;
nAvgs=floor(nTrials/nBuffers); % averaging bins
Base.NTubes=TubeList.General.NumberOfTests;
Base.NPairsTubes=...
  factorial(Base.NTubes)/(factorial(Base.NTubes-2)*factorial(2));
estTubeLengths=zeros(Base.NTubes,1);
ilag=zeros(Base.NTubes,1);
pClk=zeros(Base.NTubes,Base.nnE); % hold equiv click measured pressure waveforms (mPa)
M=-xpr.Parameters.M; % for chirp to click use minus sign (from Keefe, 2020b)
nCircEqClk=1699; % H stm important, discard initial small peaks that appear periodic
pcoeff=M*pi/Base.nd2^2;
rngff=2:Base.nd2;
theta=zeros(1,Base.nnE); % hold chirp phase
thetad2=pcoeff*(rngff-1).*(rngff-1); % sweep low to high frequencies
theta(rngff)=thetad2;
theta(Base.nd2+1)=pcoeff*Base.nd2^2;
theta((Base.nd2+2):Base.nnE)=-fliplr(thetad2); % complex conjugate to complete chirp phase
H=exp(1i*theta);
% Set up alpha for RF model
% Define tube area, radius and alpha1 for RF model
Base.TubeArea=pi/4*TubeList.General.TubeDiameter^2; % cm^2
Base.TubeRadius=sqrt(Base.TubeArea/pi);
zspecific=Base.airVec.rho*Base.airVec.c; % specific impedance of air
Base.TubeZc=zspecific/Base.TubeArea;
lv=Base.airVec.eta/zspecific; % eta/(rho*c)
lt=lv/Base.airVec.nu^2;%lt=kappa/(rho*c*Cp),Prandtl % nu=sqrt(eta*Cp/kappa)
Base.alpha1=(sqrt(lv)+Base.airVec.gamma1*sqrt(lt))/sqrt(Base.airVec.c);
chi0=1.0; % ideal viscothermal loss model chi0==1
pNorm=zeros(Base.NTubes,Base.mMat); % tube data of normalized equiv clicks
data=zeros(Base.NTubes,Base.nnE,nBuffers);
datam=zeros(Base.NTubes,Base.nnE);
dataSE=zeros(Base.NTubes,Base.nnE);
MRT=cell(1,Base.NTubes);
swVerbose=1;
nOnset=24;
hOnset=hanning(2*nOnset);
rngOnset=1:nOnset;
hOnset=hOnset(rngOnset)'; % onset windown of early click stim response

ifig=1;
for mm=[1,Base.NTubes] % loop to calculate pNorm and model RF of tube 4
  td=TubeData{mm}; % calculation is for all tubes as explained in 2020b paper
  if td.trialNum<nTrials
    disp('Bad data');
    keyboard
  end
  dataAll=td.allrecdata(:,1:nTrials);
  dataraw=reshape(dataAll,Base.nnE,nBuffers,nAvgs);
  datarawm=mean(dataraw,3);
  data(mm,:,:)=datarawm;
  cpClk=ifft(transpose(H).*fft(dataAll)); % inverse chirp filter
  pClkraw=transpose(real(cpClk));
  irng=ArtifactRejectWaveformPost(pClkraw);
  pClkallNoShift=pClkraw(irng,:); % for freq domain SNR
  pClkall=circshift(pClkallNoShift,-nCircEqClk,2); % remove leading 0's in pClkraw
  pClkkk=mean(pClkall); % mean measurement of equiv clk ADC voltage response
  % BP filter
  pClkkkPad=[zeros(1,1024),pClkkk]; % unnecessary to zero pad at end
  % Base.hhBP: impulse response function of microphone sensitivity
  temp=conv(pClkkkPad,Base.hhBP,'same');  
  temp2=circshift(temp,-1024);
  pClk(mm,:)=temp2(1:Base.nnE); % BP filtered pressure data
  pClk(mm,rngOnset)=hOnset.*pClk(mm,rngOnset); % onset ramp of 0.5 ms
  datam(mm,:)=mean(dataAll(:,irng),2);
  dataSE(mm,:)=std(dataAll(:,irng),[],2)/sqrt(sum(irng));
  ilag(mm)=PeriodicityDetector_RF(Base,pClk(mm,:)); % estimate round-trip delay in samples
  estTubeLengths(mm)=Base.airVec.c/(2*Base.fs)*ilag(mm); % in cm
  [MRT{mm},ifig]=R_ClosedTube_Time(mm,Base,estTubeLengths(mm),chi0,...
    swVerbose,ifig); % calculate model RF
  if mm==1 % need mm==1 for long-tube data to normalize in pNorm
    [Base.pMax1,Base.ipMax1]=max(abs(pClk(1,Base.ixMat)));
  end
  pNorm(mm,:)=transpose(pClk(mm,Base.ixMat))/Base.pMax1; % max ampl of 1
end % mm loop

coDefault=colororder;
hf99=figure(99);
set(hf99,'Position',[100,100,760,320]);mm=4;
subplot(1,2,1),
MRTmm=MRT{mm}
Lmm=estTubeLengths(mm);
ylim=[-0.5,16];
CalcNewModelRF_Share(Base,Lmm,MRTmm.A,ylim,coDefault);
%set(gca,'XLim',[5.2,5.45],'XTick',5.2:0.05:5.45,'YLim',[-0.5,16]);
set(gca,'XLim',[5.2,5.4],'XTick',5.2:0.05:5.4,'YLim',ylim);
hold on;
hp=plot(1e3*Base.t,1e-3*MRTmm.rcmDec,'-.*');
hp.Color=coDefault(3,:);
plot(1e3*MRTmm.tm,1e-3*MRTmm.rcm,'-k','DisplayName','Exact');
hl=legend('LW','Approximate','\tau','Decimation','Exact',...
  'Location','NorthEast','box','on','AutoUpdate','off');
plot(1e3*MRTmm.tmax,1e-3*MRTmm.rcmax,'kx');

subplot(1,2,2),
L2=2; % 2 cm tube
[MRT2,ifig]=R_ClosedTube_Time(1,Base,L2,chi0,swVerbose,ifig);
ylim=[-40,100];
CalcNewModelRF_Share(Base,L2,MRT2.A,ylim,coDefault);
set(gca,'XLim',[0.07,0.3],'YLim',ylim);
hold on;
plot(1e3*(2*MRT2.Asq/3+2*L2/Base.airVec.c),20*log10(1e-3*MRT2.rcmax),...
  'kx','DisplayName','Exact');
legend('LW','Approximate','\tau','Exact','Location','NorthEast','box','on');

% Increase font sizes by factor of 1.2
h=findall(gcf, '-property', 'fontsize');
hFontSize = cell2mat(get(h,'FontSize'));
newFontSize = hFontSize * 1.2;
set(h,{'FontSize'}, num2cell(newFontSize))
MRTmm


