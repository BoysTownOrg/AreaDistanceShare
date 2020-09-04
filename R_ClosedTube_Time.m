function [MRT,ifig]=R_ClosedTube_Time(kk,Base,Len,chi,swVerbose,ifig)
%R_ClosedTube_Time, 1/18/19. Return reflectance tube model in time
% domain.
%R_ClosedTube_RFrev, 1/7/19. Revised to use Base as input
%R_ClosedTube, Calculate reflection function of kk-th closed cylindrical 
%  tube of area (area) and length Len) with viscothermal coeff chi.
% 10/29/18. Len and chi vary in each iteration.
tau=2*Len/Base.airVec.c;
A=chi*Base.alpha1*Len/Base.TubeRadius;
Asq=A*A;
% MRT is Model reflection function struct
if isempty(ifig)
  swExact=0;
else
  swExact=1;
end
[MRT.rc0,MRT.rc0Integral,MRT.rcfun]=CalcReflectionFcn(Base.t,tau,A,Asq,...
  Base.T,swExact);
%MRT.mult=10; % Use multi-state interpolator if mult>13 
MRT.mult=15; % 3 is more accurate peak than 10 for resampling, but more undershoot <2L/c
tm0=interp([0,Base.t],MRT.mult);
MRT.tm=tm0(2:end);
[MRT.rcm,MRT.rcmIntegral]=CalcReflectionFcn(MRT.tm,tau,A,Asq,...
  Base.T/MRT.mult,0);
MRT.tmax=2/3*Asq+tau;
MRT.rcmax=(3/2)^(3/2)*exp(-3/2)/sqrt(pi)/Asq;
MRT.tau=tau;
MRT.A=A; % coefficient A not matrix A!
MRT.Asq=Asq;
% decimate
rcmM=[0,MRT.rcm;]; % add in initial 0 at time 0
MRT.rcmDec=resample(rcmM,tm0,Base.fs,1,MRT.mult,Base.hhLP,'spline');
if length(MRT.rc0)==length(MRT.rcmDec)-1 
  MRT.rcmDec=MRT.rcmDec(2:end); % remove any spurious value at 0
end
if swVerbose && swExact
  MRT.rcmDecIntegral=sum(MRT.rcmDec)*Base.T;
  MRT.rcmIntegral=sum(MRT.rcm)*Base.T/MRT.mult;
  disp(' ');
  disp(['Tube ',int2str(kk),':']);
  disp(['Numerical integral of RF over measurement duration: ',...
    num2str(MRT.rcfun,'%7.4f')]);
  disp(['Sum of original audio-rate RF: ',...
    num2str(MRT.rc0Integral,'%7.4f')]);
  disp(['Sum of RF upsampled by ',int2str(MRT.mult),': ',...
    num2str(MRT.rcmIntegral,'%7.4f')]);
  disp(['Sum of RF upsampled/LPfiltered/downsampled to audio rate: ',...
    num2str(MRT.rcmDecIntegral,'%7.4f')]);
  if kk==4
    temp=1e-3*MRT.rcmax;
    disp(['Max cts-time RF: ',num2str(temp,'%7.4f'),' 1/ms']);
    disp(['Max RF upsampled by ',int2str(MRT.mult),': ',...
      num2str(1e-3*max(MRT.rcm),'%7.4f'),' 1/ms']);
    disp(['Max audio-rate RF: ',num2str(1e-3*max(MRT.rc0),'%7.4f'),' 1/ms']);
    disp(['Max RF upsampled/LPfiltered/downsampled to audio rate: ',...
      num2str(1e-3*max(MRT.rcmDec),'%7.4f'),' 1/ms']);
    disp(['Min cts-time RF: ',num2str(0,'%7.4f'),' 1/ms']);
    disp(['Min RF upsampled by ',int2str(MRT.mult),': ',...
      num2str(1e-3*min(MRT.rcm),'%7.4f'),' 1/ms']);
    disp(['Min audio-rate RF: ',num2str(1e-3*min(MRT.rc0),'%7.4f'),' 1/ms']);
    disp(['Min RF upsampled/LPfiltered/downsampled to audio rate: ',...
      num2str(1e-3*min(MRT.rcmDec),'%7.4f'),' 1/ms']);
  end
end
end

function [rc,rcIntegral,rcfun]=CalcReflectionFcn(t,tau,A,Asq,T,swExact)
% Calculate model reflection function of closed tube
rc=zeros(size(t));
ix=find(t-tau>0);
ttauinv=1./(t(ix)-tau);
rc(ix)=(A/sqrt(pi)*ttauinv.^(1.5)).*exp(-Asq.*ttauinv);
if swExact
  rcIntegral=T*sum(rc);
  tend=t(end);
  fun=@(x) (A/sqrt(pi)./(x-tau).^(1.5)).*exp(-Asq./(x-tau));
  rcfun=integral(fun,tau,tend,'RelTol',1e-8,'AbsTol',1e-13);
else
  rcIntegral=0;
  rcfun=0; % not calculated
end
end
