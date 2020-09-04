function ilag=PeriodicityDetector_RF(Base,pt)
% PeriodicityDetectorRF for reflection function modified from
% PeriodicityDetector in Lark50.
% Periodicity detector based on autocorrelation function using the Sondhi
% center clipped model.  Reference"Digital Processing of Speech Signals",
% Rabiner&Schafer (Prentice Hall, Englewood Cliffs, 1978), pages 150-156.

c=Base.airVec.c; % cm/s
sRate=Base.fs;
lenmin=80; % min tube length in cm
lenmax=200; % min tube length in cm
minlag=2*lenmin*sRate/c;
maxlag=2*lenmax*sRate/c;
[ct,Nmin,Nmax]=CenterClip(pt);
Nr=min(Nmax,round(1.1*maxlag)); % need at least maxlag delays
rn=zeros(1,Nr); % holds autocorrelation function with zero delay at kk=1
for kk=1:Nr
  for mm=Nmin:Nmax+1-kk
    if ct(mm)~=0 && ct(mm+kk-1)~=0 % calculate only for non-clipped samples
      rn(kk)=rn(kk)+ct(mm)*ct(mm+kk-1);
    end
  end
end

%dimax=4; % dimax>=2 to omit values near origin at sample 1
dimax=10; % dimax>=2 to omit values near origin at sample 1
[~,imax]=max(rn(dimax:Nr)); % don't look at initial samples near 0 delay
nz=3; % number of consecutive zeros to look for
if imax<=round(minlag)
  jj=1;
  while (jj<Nr-nz+1) && max(abs(ct(jj:jj+nz)))==0
    jj=jj+1;
  end
  if (jj>=Nr-nz+1)
    ilag=0;
    disp('No solution found');
    keyboard
  else
    [~,imax2]=max(rn(jj:Nr));
    ilag=imax2+jj-2; %-1 to balance jj==1, and another -1 to remove 0 delay
  end
else
  ilag=imax+dimax-2; 
end

% make quadratic fit to peak of ACF to refine estimate of peak delay
% yy=a*xx.^2+b*xx+c
xx=[-1 0 1];
yy=rn((ilag):(ilag+2)); %rn(1)==0 ms delay, so the peak is offset by 1 sample
% this is the explicit quadratic solution for equally-spaced samples xx
a=0.5*yy(1)-yy(2)+0.5*yy(3);
b=-0.5*(yy(1)-yy(3));
xpeak=-0.5*b/a; % lies between +/-0.5
ilag=ilag+xpeak; % ilag is no longer an integer

swPlot=0;
if swPlot==1
  tt=1:Nmax;
  ttr=0:(Nr-1);
  hf1=figure;
  subplot(2,1,1)
  plot(tt,pt(1:Nmax));
  ylabel('Raw waveform');
  if hf1==1
    ss='First';
  else
    ss='Second';
  end
  title([ss, ' tube waveform'])
  subplot(2,1,2)
  plot(tt,ct(1:Nmax));
  xlabel('Time (samples)');
  ylabel('Clipped waveform');
  figure;
  yabs=max(abs(rn));
  plot(ttr,rn,'-b',ilag*[1 1],yabs*[-1 1],'--k');
  xlabel('Delay (samples)');
  ylabel('Autocorrelation function');
  title([ss, ' tube ACF'])
  % check Matlab solution, no need to implement this except diagnostics
  p=polyfit(xx,yy,2);
  x=-1:0.02:1;
  y=polyval(p,x);
  c2=yy(2); % not needed except in plot
  figure;
  plot(x,y,'-b',xx,yy,'bd',xpeak,a*xpeak^2+b*xpeak+c2,'ko');
  title([ss, ' peak fitting'])
  xlabel('Delay re: peak (samples)');
  ylabel('ACF near peak');
end
end

function [ct,Nmin,Nmax]=CenterClip(pt)
% Perform center clipped waveform using Sondhi method
maxamp=max(abs(pt));
clip=0.01*maxamp;
ct=zeros(size(pt));
ip=pt>clip;
im=pt<-clip;
ct(ip==1)=pt(ip==1)-clip;
ct(im==1)=pt(im==1)+clip;
temp=find(ip==1);
ipos=max(temp);
ipos0=min(temp);
temp=find(im==1);
ineg=max(temp);
ineg0=min(temp);
Nmax=max(ipos,ineg);
Nmin=min(ipos0,ineg0);
end
