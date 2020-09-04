classdef BaseRFShare<handle
  % BaseRFShare is base class for AreaDistanceShare applications 
  % 8/14/20. Code adapted from BaseRF
  
  properties
    gprm
    uprm
    xpr
    xprData
  end
  
  properties(GetAccess=public, SetAccess = private)
    fs=48000;
    fskHz
    freqs
    freqskHz
    T
    Tms
    AnalysisFreqLow=176.7767; % H stm, 1/2th octave below 0.25 kHz
    AnalysisFreqHi=14254.4; % H stm, 1/6th octave below 16 kHz
    hhBP % BP filter impulse response to filter measured pressures
    nhhBP % length of hhBP
    hhLP % LP filter impulse response to resample RFs
    nhhLP % length of hhLP
    nPmultiple=3.75 % mMat1 for single tube is nP=nPmultiple*1024
    nP
    mMat % # of samples/rows in A matrix eq./tube
    md2
    ixMat
    t
    ff % full-range frequencies based on fs/mMat
    ixfAnalysis % frequency range in analysis band
    fAnalysis % frequencies (Hz) in analysis band
    ixfAnalysisUpper % frequency range in upper half DFT > 24 kHz & <48 kHz
    logf % log2 of fAnalysis
    NfAnalysis
  end
  
  properties(GetAccess=public, SetAccess = public)
    % for buffer handling and analyses
    nnE
    nd2
    NTubes % used in analyses
    NPairsTubes
    altitude
    temperature
    alpha1
    TubeRadius
    TubeArea
    TubeZc % characteristic impedance of tube, real in lossfree limit
    pMax1; % max incident clk pressure (Tube 1) used to normalize all pNorm 
    ipMax1; % sample number of max pressure in Tube 1
  end
  
  properties(Dependent)
    airVec
  end
  
  %% public methods
  methods
    %% Constructor
    function self = BaseRFShare(swRealTime)
    % swRealTime is unused
    end   
        
    function airvec=get.airVec(self)
      % AirParams2, 11/30/09. Output CGS thermodynamic constants of air.
      % [Based on AH Benade, J. Acoust. Soc. Am. (1968) with exponential
      %     atmosphere model added]
      %   Atmospheric pressure and bulk modulus of air use altitude 
      %   correction (SI units) on rho and Pressure
      % class Inputs: Temperature (deg C) and altitude (m), 
      %    rounding altitude to nearest 100 m suffices.]
      dT=self.temperature-26.85;  %Temperature is in degrees C
      % All constants in CGS units in the main function
      airvec.c=34723*(1.0+0.00166*dT); % independent of altitude
      rho0=0.0011769*(1.0-0.00335*dT); % at sea level, altitude=0
      % Altitude correction for pressure and density using
      % exponential atmosphere model for dry air on earth under isothermal
      % conditions.  OK because test room has const. temperature as
      % altitude varies.
      RgasconstantMKS=8.31432; % universal gas constant, J/(mole*K)
      g0MKS=9.80665; % acceleration of gravity, m/s .
      Mm=0.0289644; % molar mass of air, kg/mole, 
      %   about 29 g for 1 mole of 1 cubic meter of air
      % AltitudeFactor is dimensionless, and equals 1 for altitude=0 m.
      AltitudeFactor=...
        exp(-g0MKS*Mm*self.altitude/(RgasconstantMKS*(self.temperature...
        +273.15)));
      airvec.rho=rho0*AltitudeFactor; % adjust rho for altitude
      airvec.eta=1.846*10^(-4)*(1.0+0.0025*dT); % independent of altitude
      airvec.nu=0.841*(1.0-0.0002*dT); % independent of altitude,nu=sqrt(eta*C_p/kappa)
      gamma=1.4017*(1-0.00002*dT); % independent of altitude
      airvec.gamma1=gamma-1;
      airvec.BulkModulus=airvec.rho*airvec.c^2; % bulk modulus of air
      airvec.Pressure=airvec.BulkModulus/gamma; % air pressure at this temperature and altitude
    end
    
    function value=get.nd2(self)
      value=fix(self.nnE/2);
    end
    
    function value=get.fskHz(self)
      value=1e-3*self.fs;
    end
    
    function myfreqs=get.freqs(self)
      myfreqs=(0:self.nd2)*self.fs/self.nn;
    end
    
    function myfreqskHz=get.freqskHz(self)
      myfreqskHz=(0:self.nd2)*self.fskHz/self.nn;
    end
    
    function value=get.T(self)
      value=1/self.fs;
    end
    
    function value=get.Tms(self)
      value=1/self.fskHz;
    end
    
    function value=get.hhBP(self) % FIR Kaiser bandpass filter
      if self.AnalysisFreqHi>15000
        fmax=self.AnalysisFreqHi+1000;
      else
        fmax=15000;
      end
      fbw=[10,self.AnalysisFreqLow,self.AnalysisFreqHi,fmax];
      % Design BP filter to convolve with band-limited measured pressures. 
      a=[0 1 0];
      r=[0.1, 0.01, 0.001];  % ripple in stop, pass, stop
      [nk,Wn,beta,ftype] = kaiserord(fbw,a,r,self.fs);
      nk = nk + rem(nk,2); %even nk, odd length(hh)
      value=fir1(nk,Wn,ftype,kaiser(nk+1,beta),'noscale');
    end
    
    function value=get.nhhBP(self)
      value=length(self.hhBP);
    end

    function value=get.hhLP(self) % FIR Kaiser lowpass filter
      % Much shorter LP filter than BP filter to resample model RFs
      if self.AnalysisFreqHi>17000
        fmax=self.AnalysisFreqHi+1000;
      else
        fmax=17000;
      end
      fbwLP=[self.AnalysisFreqHi,fmax]; % FIR Kaiser bandpass filter
      aLP=[1 0];
      rLP=[0.01, 0.01];  % ripple in stop, pass, stop
      [nkLP,WnLP,betaLP,ftypeLP] = kaiserord(fbwLP,aLP,rLP,self.fs);
      nkLP = nkLP + rem(nkLP,2); %even nk, odd length(hh)
      value=fir1(nkLP,WnLP,ftypeLP,kaiser(nkLP+1,betaLP),'noscale');
    end

    function value=get.nhhLP(self)
      value=length(self.hhLP);
    end
    
    function value=get.nP(self)
      value=self.nPmultiple*1024;
    end
    
    function value=get.mMat(self)
      value=min(self.nnE,self.nP);
    end
    
    function value=get.md2(self)
      value=fix(0.5*self.mMat);
    end
    
    function value=get.ff(self)
      value=(0:self.md2)*self.fs/self.mMat;
    end
    
    function value=get.ixfAnalysis(self)
      value=find(self.ff>=self.AnalysisFreqLow & self.ff<=self.AnalysisFreqHi);
    end
    
    function value=get.ixfAnalysisUpper(self) % indices for upper half of DFT
      value=self.mMat+2-self.ixfAnalysis;
    end
    
    function value=get.fAnalysis(self)
      value=self.ff(self.ixfAnalysis);
    end
    
    function value=get.logf(self)
      value=log2(1e-3*self.fAnalysis); % for plotting re:0 at 1 kHz
    end
    
    function value=get.NfAnalysis(self)
      value=length(self.ixfAnalysis);
    end
    
    function value=get.ixMat(self)
      value=1:self.mMat;
    end
    
    function value=get.t(self)
      value=self.T*(self.ixMat);
    end
    
  end
end


