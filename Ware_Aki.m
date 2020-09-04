function [r,area]=Ware_Aki(T,rf,area0,Nt)
% Ware_Aki, 4/3/19.
% Ware_Aki (1969) algorithm to calculate area-distance function for unknown 
%  waveguide based on Sharp dissertation
% Inputs:
%   rf, reflection function of unknown waveguide at each increment of 
%     sample period T
%   area0, area at entryway
%   Nt, number of samples to calculate area <= length(rf)
% Outputs:
%   area, area distance function as function of air-layer location z
%   r, reflectance as function of z, delta z = c*T/2 for sample period T

%radius: starting radius (m)
%finlen: final length to reconstruct up to
%c: speed of sound (m/s)
%samp_num: length of iir
%time: time in ms
% input RF array sampled at times t=1:samp_num
if Nt<2
  disp('Choose Nt>1');
  return;
elseif Nt>length(rf)
  disp('RF too short to estimate Nt samples of area');
  return
end
rfDimensionless=T*rf;
% Sharp uses times from from 0 to samp_num, use 1:samp_num here
% Init variables
product=1;
F=zeros(1,Nt);
G=zeros(1,Nt);
r=zeros(1,Nt); % reflectance of successive air layers
F(1)=1;
area=zeros(1,Nt);
areaCurrent=area0;
for ii=1:Nt
  if ii>1
    product=product*(1-r(ii-1))*(1+r(ii-1));
  end
  numerator=0;
  for jj=1:ii
    numerator=numerator+F(jj)*rfDimensionless(ii-jj+1);
  end
  r(ii)=numerator/product;
  areaCurrent=areaCurrent*(1-r(ii))/(1+r(ii));
  area(ii)=areaCurrent;
  if areaCurrent<0
    disp(['ii=',int2str(ii),': area cannot be negative']);
    return
  end
  jj=ii;
  while jj>=2 % recursive update of F and G
    F(jj)=F(jj)+r(ii)*G(jj-1);
    G(jj)=r(ii)*F(jj)+G(jj-1);
    jj=jj-1;
  end
  F(1)=1;
  G(1)=r(ii);
end
r=r/T; % convert to units of (1/time) based on units of T
