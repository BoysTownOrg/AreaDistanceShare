function [r,area]=LayerPeeling_Cylinder_Lossless(T,rf,area0,Nt)
% LayerPeeling_Cylinder_Lossless, 6/7/19.
% Based on AM Bruckstein and T Kailath, Inverse scattering for discrete
%   transmission-line models, SIAM Review 29, 359-389 (1987);
% as modified for forward/reverse pressure by
% N Amir, U Shimony and G Rosenhouse, A discrete model for tubular
% acoustic systems with varying cross section--direct and inverse problems,
%  Acustica 81, 450-462 (1995);
% as corrected and porting C code from
% DB Sharp, Ph.D dissertation, "Acoustic pulse reflectometry for the
% measurement of musical wind instruments ," U. Edinburgh (1996).

% Inputs:
%   T, sample period
%   rf, reflection function of unknown waveguide at each increment of
%     sample period T
%   area0, area at entryway
%   Nt, even number of samples to calculate area <= length(rf)
% Outputs:
%   area, area distance function as function of air-layer location z
%   r, reflectance as function of z, delta z = c*T/2 for sample period T
if Nt<2
  disp('Choose Nt>1');
  return;
elseif Nt>length(rf)
  disp('RF too short to estimate Nt samples of area');
  return
end

input_real=zeros(1,Nt); % values of input-real time shift in main loop
input_real(1)=1;
reflec_real=rf*T; % values of rf time shift in reflec_real in main loop
r=zeros(1,Nt); % reflectance of successive air layers
area=zeros(1,Nt); % area of tube for each air layer
areaCurrent=area0;
for ii=1:Nt % main loop peeling off each layer
  r1=reflec_real(1)/input_real(1);
  r(ii)=r1;
  areaCurrent=areaCurrent*(1-r1)/(1+r1); % recursive update of area
  area(ii)=areaCurrent;
  if areaCurrent<0
    disp(['ii=',int2str(ii),': area cannot be negative']);
  end
  denom1r1=1/(1-r1);
  input_real_old=input_real;
  input_real=(input_real_old-r1*reflec_real)*denom1r1; %update forward pres
  reflec_real=(reflec_real-r1*input_real_old)*denom1r1;%update reverse pres
  reflec_real=circshift(reflec_real,-1); % delay reverse wave by 1 sample
  reflec_real(Nt)=0;
end
r=r/T; % convert to units of (1/time) based on units of T
