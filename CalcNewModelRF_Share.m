function CalcNewModelRF_Share(Base,L,A,ylim,coDefault)
% CalcNewModelRF_Share, 8/17/20.  Calculate cylindrical tube model RF in
% time domain based on LW approach. Exact and approximation RF's.
tauL=2*L/Base.airVec.c;
T=Base.T;
t0=Base.t; % starts at t=T,2T, etc.
ixTau=t0>tauL;
ixTauTD2=t0>tauL+T/2;

temp1=zeros(1,Base.mMat);
temp2=zeros(1,Base.mMat);
temp1(ixTau)=1-erf(A./sqrt(t0(ixTau)-tauL+0.5*T));
temp2(ixTauTD2)=1-erf(A./sqrt(t0(ixTauTD2)-tauL-0.5*T));
rNew=temp1-temp2;

tTaylor=t0(ixTauTD2)-tauL;
Asq=A*A;
rTaylor1=1/(sqrt(pi))*A.*tTaylor.^(-1.5).*exp(-Asq./tTaylor); % 3/1/20
if L<50 % short tubes log plot
  ep=1e-14; % small value as lower limit to define dB
  ixx=rNew<ep;
  rNew(ixx)=repmat(ep,1,sum(ixx));
  LrNew=20*log10(1e-3*rNew/T);
  ixx=rTaylor1<ep;
  rTaylor1(ixx)=repmat(ep,1,sum(ixx));
  LrTaylor1=20*log10(1e-3*rTaylor1);
  plot(1e3*t0,LrNew,'-o',1e3*(tTaylor+tauL),LrTaylor1,'--');
  ylabel('RF level (dB),  0 dB re: 1/ms');
  text(0.14,LrNew(6),'n=6 at time {\itnT}',...
    'HorizontalAlignment','left','VerticalAlignment','middle');
else
  plot(1e3*t0,1e-3*rNew/T,'-o',1e3*(tTaylor+tauL),1e-3*rTaylor1,'--');
  ylabel('RF (1/ms)');
end
hold on;
hp=plot(1e3*[tauL,tauL],ylim,':');
hp.Color=coDefault(4,:);
hp.LineWidth=0.7;
xlabel('Time (ms)');
title(['Length = ',sprintf('%4.1f',L),' cm']);

sum(rNew) % sum of RF
end

