function LossyPropagation(Base,MnearMax,sL,gl,gr,areaCurrent)
% 4/23/20. Uses 0th order Taylor series for q0, 1st order for qk.
swLossy=1;
T=Base.T;
Tms=Base.Tms;
s=Base.airVec;
ifig=1;
aT=sqrt(T/2)/2*sqrt(s.eta/s.rho)*(1+s.gamma1/s.nu); % in cm
% constant area here, but not in general loop
radiusCurrent=sqrt(areaCurrent/pi);
x0=aT/radiusCurrent; % used in PropagateP,M
% Evaluate accuracy over Naf identical segments
MnearMaxM1=MnearMax-1;
nOffset=MnearMaxM1/2;
% grF=zeros(MnearMax,Naf); % forward right normalizedpressure, index m=1 at input m=0, etc.
% grF(1,1)=1;
nT=(0:MnearMaxM1)+nOffset;
tn=[0,nOffset-1,nT]*Tms;
%swPositive=1; %=0 path is no good because stimulus is at m=0, not MnearMax
% Will have to test PropagateM within area estimation loop later
%if swPositive
  glP=gl; % gl and gr are all zeros, so no big deal
  grP=gr;
  grP00=1; % forward impulse
  glP(1,1)=grP00; % no discontinuity at m=1.
  for m=1:MnearMax % loop for cylindrical tube
    %for m=1:5 % loop for cylindrical tube
    nrng1=1:m;
    grP(m,nrng1)=PropagateP(m,glP(m,nrng1),x0,swLossy); % grP is 1/2 T later at same n than glP
    glP(m+1,nrng1)=grP(m,nrng1);
  end
  glP(1:9,1:10)
  grP(1:9,1:10)
  glF=glP; % Here work with glP, not grP
% else % solve using calls to PropagateM in same tube
%   glM=gl;
%   grM=gr;
%   glM00=1; % incident reverse impulse, actually not 0,0 but opposite end, and l/r would be interchanged too
%   glM(1,1)=glM00; % no discontinuity at m=1.
%   for m=1:MnearMax % loop for cylindrical tube
%     %for m=1:5 % loop for cylindrical tube
%     nrng1=1:m;
%     grM(m,nrng1)=PropagateM(m,glM(m,nrng1),x0,swLossy); % grP is 1/2 T later at same n than glP
%     glM(m+1,nrng1)=grM(m,nrng1);
%   end
%   glM(1:9,1:10)
%   grM(1:9,1:10)
%   glF=glM; % Here work with glM, not grM
% end

rn=[0,0,(glF(MnearMax,1:MnearMax))/Tms];
plot(tn,rn,'-o','DisplayName','LW');
hold on;
disp(['Current integral of RF is: ',sprintf('%7.5f',Tms*sum(rn))]);
title(sL);
set(gca,'XLim',[max(nOffset*Tms-0.15,0),nOffset*Tms+0.2],'YLim',[-1,51]);
% Construct tube transmission model
Base.TubeRadius=sqrt(Base.TubeArea/pi);
zspecific=Base.airVec.rho*Base.airVec.c; % specific impedance of air
Base.TubeZc=zspecific/Base.TubeArea;
lv=Base.airVec.eta/zspecific; % eta/(rho*c)
lt=lv/Base.airVec.nu^2;%lt=kappa/(rho*c*Cp),Prandtl % nu=sqrt(eta*Cp/kappa)
Base.alpha1=(sqrt(lv)+Base.airVec.gamma1*sqrt(lt))/sqrt(Base.airVec.c);
Len=MnearMaxM1*(Base.airVec.c*T/2);
[MRT,~]=R_ClosedTube_Time(1,Base,Len/2,1,1,ifig);
hold on;
nT0=Tms*(1:MnearMax); % slide to align
hp2=plot(nT0,1e-3*MRT.rcmDec(1:MnearMax),'--x','DisplayName','Decimation');

ntemp=MRT.mult*MnearMax;
nT0M=(0:(ntemp-1))*Tms/MRT.mult;
plot(nT0M,1e-3*MRT.rcm(1:ntemp),'DisplayName','Nearly Exact');
xlabel('Time (ms)');
ylabel('RF (1/ms)');
hl=legend;
set(hl,'Box','on');

  function [grPm,x0]=PropagateP(m,glPm,x0,swLossy)
    % Calculate convolution product out to n=m
    % This retains highest power of x0 to be x0^2
    if swLossy==0
      grPm=glPm; % time delay handled in calling code
    else
      x0pi=x0/sqrt(pi);
      N=length(glPm);
      grPm=zeros(1,N);
      A0=1-2*x0pi;
      grPm(1)=A0*glPm(1);
      % Index of time is grPm and glPm are each adjusted so that index 1 of
      % MATLAB array denotes time on main diagonal of space-time diagram.
      % Thus, the half-sample delay is not apparent in these array pairs.
      if m>1
        rng1=2:m;
        sum2=zeros(1,m-1); % 1st term always zero
        for n=rng1
          sumk=0; % zero for each n
          for k=1:(n-1) % always n-k>0 within this loop
            r2km1=sqrt(2*k-1);
            r2kp1=sqrt(2*k+1);
            Itildek0=2*(1/r2km1-1/r2kp1);
            Itildek1=r2kp1-r2km1;
            temp1=glPm(n-k)*((1+k)*Itildek0-Itildek1); % * ck term
            if n-k-1>0
              temp1=temp1+...
                glPm(n-k-1)*(-k*Itildek0+Itildek1); % * dk term
            end
            sumk=sumk+temp1;
          end
          sum2(n-1)=sumk;
        end
        grPm(rng1)=A0*glPm(rng1)+x0pi*sum2;
      end
    end
  end % function

  function [grMm,x0]=PropagateM(m,glMm,x0,swLossy)
    % This retains only linear term in x0. Not yet tested as of 2/19/20 but
    % appears parallel to code in PropagateP
    if swLossy==0
      grMm=glMm; % time delay handled in calling code
    else
      x0pi=x0/sqrt(pi);
      N=length(glMm);
      grMm=zeros(1,N);
      A0=1-2*x0pi;
      grMm(1)=glMm(1)/A0;
      if m>1
        rng1=2:m;
        sum2=zeros(1,m-1); % 1st term always zero
        for n=rng1
          sumk=0; % zero for each n
          for k=1:(n-1) % always n-k>0 within this loop
            r2km1=sqrt(2*k-1);
            r2kp1=sqrt(2*k+1);
            Itildek0=2*(1/r2km1-1/r2kp1);
            Itildek1=r2kp1-r2km1;
            temp1=glMm(n-k)*((1+k)*Itildek0-Itildek1); % * ck term
            if n-k-1>0
              temp1=temp1+...
                glMm(n-k-1)*(-k*Itildek0+Itildek1); % * dk term
            end
            sumk=sumk+temp1;
          end
          sum2(n-1)=sumk;
        end
        grMm(rng1)=(glMm(rng1)-x0pi*sum2)/A0;
      end
    end
  end

end

