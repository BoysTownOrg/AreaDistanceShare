function [grP,grM,glP,glM,areaD,rD]=CalcNearTM2(Base,swLossy,Naf,hRF,...
  MnearMax)
% CalcNearTM2, 1/25/22. When called for ear-referenced quantities, the
% area at the initial sample for r1 is not Base.TubeArea but rather the
% initial area determined in the first deconvolution.
% CalcNearTM, 2/25/2020.  Calculate area-distance function using layer
% peelingfrom measured RF at probe tip (hRF) over duration of Naf samples
% for discontinuities from 1:MnearMax with or without a model of
% viscothermal loss (swLossy).

s=Base.airVec;
aT=sqrt(Base.T/2)/2*sqrt(s.eta/s.rho)*(1+s.gamma1/s.nu); % in cm
% Spatial index m goes from 1 to Mnear, so first MATLAB array index goes
%   from 1 to Mnear.
% Time index n goes from 0 to N-1 for m=1.  For all areas with m>1,
%  first MATLAB array index is at n=0 relative to the earliest time at
%   which the normalized pressure at m is non-zero, i.e., for
%   g_r(m-1,(m-1)/2)~=0. See Fig. 1 of Keefe (2020a). So times in second 
%   index are delayed relative to time zero.
glP=zeros(MnearMax,Naf); % g_l^+
glM=zeros(MnearMax,Naf); % g_l^-
grP=zeros(MnearMax,Naf); % g_r^+[m,n] at time 1/2 sample later than
%  g_l^+[m,n] so g_r^+[m,n] represents g_r^+[m,n+1/2]
grM=zeros(MnearMax,Naf); % g_r^-[m,n] at time 1/2 sample earlier than
%  g_l^-[m,n] so g_r^-[m,n] represents g_r^-[m,n-1/2]
areaD=zeros(MnearMax,1); % Area at discontinuity 1:MnearMax
radiusD=zeros(MnearMax,1); % Radius in segment at discontinuity 1:MnearMax
rD=zeros(MnearMax,1); % starts at m=1, RF at discontinuity m

% Init for discontinuity between m=0 and m=1,
r1=hRF(1); % time 0 is hRF(1)
rD(1)=r1; % RF at m==1 for entryway of ear
areaCurrent=Base.TubeArea*(1-r1)/(1+r1); % area of entry tube segment 1
areaD(1)=areaCurrent;
radiusCurrent=sqrt(areaCurrent/pi); % in cm
radiusD(1)=radiusCurrent;
x0pi=aT/(radiusCurrent*sqrt(pi)); % used in PropagateP,M
grP00=1; % forward impulse
grM00=hRF(1); % forward measured RF
denom1r1=1/(1-r1);

glP(1,1)=denom1r1*(grP00-r1*grM00);
nrng2=2:Naf;
glP(1,nrng2)=denom1r1*(-r1*hRF(nrng2));
glM(1,nrng2)=denom1r1*hRF(nrng2);
nrng1=1:(Naf-1);
Nmax=Naf;
for m=2:MnearMax % go further with MnearMax > expected largest m
  grP(m-1,nrng1)=PropagateP(glP(m-1,nrng1),Nmax-1,x0pi,swLossy);
  % grP is 1/2 T later at same n than glP
  grM(m-1,nrng1)=PropagateM(glM(m-1,nrng2),Nmax-1,x0pi,swLossy);
  % glM is 1/2 T earlier at same n than grM
  r1=grM(m-1,1)/grP(m-1,1); % same as r1 above for m=1, new for m>1
  rD(m)=r1;
  areaCurrent=areaCurrent*(1-r1)/(1+r1); % recursive update of area
  areaD(m)=areaCurrent;
  radiusCurrent=sqrt(areaCurrent/pi);
  radiusD(m)=radiusCurrent;
  x0pi=aT/(radiusCurrent*sqrt(pi)); % used in PropagateP,M
  denom1r1=1/(1-r1);
  glP(m,nrng1)=denom1r1*(grP(m-1,nrng1)-r1*grM(m-1,nrng1));
  glM(m,nrng1)=denom1r1*(-r1*grP(m-1,nrng1)+grM(m-1,nrng1));
  Nmax=Nmax-1;
  nrng1=1:(Nmax-1);
  nrng2=2:Nmax;
end % areaD in this loop agree with Sharp code for areaD and rD
%   at node m from 1 to MnearMax

  function grPm=PropagateP(glPm,N,x0pi,swLossy)
    % Calculate convolution product out to n=m
    % This retains highest power of x0 to be x0^2
    if swLossy==0
      grPm=glPm; % time delay handled in calling code
    else
      grPm=zeros(1,N);
      A0=1-2*x0pi;
      grPm(1)=A0*glPm(1);
      % Index of time is grPm and glPm are each adjusted so that index 1 of
      % MATLAB array denotes time on main diagonal of space-time diagram.
      % Thus, the half-sample delay is not apparent in these array pairs.
      rng1=2:(N-1);
      sum2=zeros(1,N-2); % 1st term always zero
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

  function grMm=PropagateM(glMm,N,x0pi,swLossy)
    % This retains only linear term in x0. Parallel to code in PropagateP
    if swLossy==0
      grMm=glMm; % time delay handled in calling code
    else
      grMm=zeros(1,N);
      A0=1-2*x0pi;
      grMm(1)=glMm(1)/A0;
      rng1=2:(N-1);
      sum2=zeros(1,N-2); % 1st term always zero
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

