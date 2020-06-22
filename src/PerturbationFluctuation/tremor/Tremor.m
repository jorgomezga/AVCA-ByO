function [rFftr, rFatr, rFTRI, rATRI]=Tremor( vSignal, iFs )

% Calculates the tremor (parameters that modulate the pitch period and
% amplitude of the voice pitch). It is like if the disturbances of the pitch period and
% the pitch amplitude are due to two modulating signals that modulate the period of
% average pitch and the average pitch amplitude respectively. In this manner the 
% periodicity and intensity of the disturbance is calculated
% The algorithm followed by Kay is used for the calculation.
%
% Input parameters:
%    vSignal:   is the voice section to analyze
%    iFs:       is the sample rate of the voice signal.
%
% Output parameters:
%   rFftr:      frequency of the low-frequency component of higher intensity that modulates
%               the pitch period.
%   rFatr:      frequency of the low-frequency component of higher intensity that modulates
%               the pitch amplitude.
%   rFTRI:      (Fo Tremor Intensity Index)
%   rATRI:      (Amplitude Tremor Intensity Index)

if nargin<2, iFs=25000; end

% Check that the vector is of type column 
if ~isvector(vSignal)
    error( 'El parámetro vSignal no es un vector columna' ); 
elseif size( vSignal, 2 ) ~= 1
    vSignal=vSignal'; 
end 

[~, ~, pitch_cont, amp_cont]=ParamPitch( vSignal, iFs, 1, length(vSignal), 1);

N=length( pitch_cont ); 

FoData=zeros(1,N);

for i=1:N
   if pitch_cont(i)~=0
      FoData(i)=1/pitch_cont(i);
   else
      % If the calculated pitch is zero for a sustained (incorrect) vowel,
      % takes either the value of the last non-zero pitch, or interpolate
      FoData(i)=FoData(i-1);
   end
end

AoData=amp_cont;

%==================================================================================
%========================== Preprocessing ===========================
%==================================================================================

% Low-pass filtering at 30 Hz (we use order 22 hamming window)
B=fir1( 22, 30/(iFs/2) );

FoDataFilt=filter(B,1,FoData);
AoDataFilt=filter(B,1,AoData);

% Subsampling at 400 Hz (we do it with the implementation in C in mind); previous
% filtering is already used to avoid aliasing:

% The decimation factor is
D=round(iFs/400);

% The resulting sampling frequency due to rounding is somewhat different
% to 400 Hz
iFsd=iFs/D;

% The number of samples of the subsampled samples will be
Nd=fix(N/D);

FoDataSubm= zeros(1, Nd);
AoDataSubm= zeros(1, Nd);
for i=1:Nd
   FoDataSubm(i)=FoDataFilt(i*D);
   AoDataSubm(i)=AoDataFilt(i*D);
end

%==================================================================================
%=================== Processing of all the parameters =======================
%==================================================================================

% The total energy of the resulting signal is
EFo=sum(FoDataSubm.^2);
EAo=sum(AoDataSubm.^2);

% Quitamos la componente DC. Se calculan así las señales que modulan al periodo
% de pitch y a la amplitud de los picos
% We remove the DC component. Thus the signals that modulate the  
% pitch period and the amplitude of the peaks
moduladora_Fo=FoDataSubm-mean(FoDataSubm);
moduladora_Ao=AoDataSubm-mean(AoDataSubm);

%EFo=sum(moduladora_Fo.^2);
%EAo=sum(moduladora_Ao.^2);

% Calculation of the autocorrelation of the resulting residual signals (modulating signals)
RFo=xcorr(moduladora_Fo);
RFo=RFo(Nd:2*Nd-1);

RAo=xcorr(moduladora_Ao);
RAo=RAo(Nd:2*Nd-1);

% Division by total energy and conversion to%
RFo=100*(RFo/EFo);
RAo=100*(RAo/EAo);

% We will look for the tremor in the frequency range 2 Hz to 10 Hz
fmin=2;
fmax=10;

% In samples
kmax=round(iFsd/fmin);
kmin=round(iFsd/fmax);

%==================================================================================
%==========================   Fftr y FTRI  =============================
%==================================================================================
[m,pos]=max(RFo(kmin:kmax));

PosiblePeriodo=pos+(kmin-1);
% We look at harmonic frequencies in case we have reached the second peak (or higher)
% The search margin is 20% of the period found
m2=m;
fin=0;
while m2>0.6*m && fin==0
   
   mueBusq=round(PosiblePeriodo/2);
   margen=round(0.2*PosiblePeriodo);
   mueIni=max(mueBusq-margen,kmin);
   mueFin=mueBusq+margen;
   
   if mueIni>kmin
      [m2,pos]=max(RFo(mueIni:mueFin));
      Pico=pos+(mueIni-1);
      if RFo(Pico)>mean(RFo(Pico-(kmin-1):Pico-1))
          % In case the first descent of the autocorrelation is reached
         PosiblePeriodo=Pico;
         m=m2;
      else
         fin=1;
      end
   else
      fin=1;
   end
   
end

Periodo=PosiblePeriodo-1;

rFftr=iFsd/Periodo; % Fftr in Hz
rFTRI=abs(m); % in %

%==================================================================================
%==========================   Fatr y ATRI  =============================
%==================================================================================
[m,pos]=max(RAo(kmin:kmax));

PosiblePeriodo=pos+(kmin-1);
% We look at harmonic frequencies in case we have reached the second peak (or higher)
% The search margin is 20% of the period found
m2=m;
fin=0;
while m2>0.6*m && fin==0
   
   mueBusq=round(PosiblePeriodo/2);
   margen=round(0.2*PosiblePeriodo);
   mueIni=mueBusq-margen;
   mueFin=mueBusq+margen;
   
   if mueIni>kmin
      [m2,pos]=max(RAo(mueIni:mueFin));
      Pico=pos+(mueIni-1);   
      if RAo(Pico)>mean(RAo(Pico-(kmin-1):Pico-1))
          % In case the first descent of the autocorrelation is reached
         PosiblePeriodo=Pico;
         m=m2;
      else
         fin=1;
      end
   else
      fin=1;
   end
   
end

Periodo=PosiblePeriodo-1;

rFatr=iFsd/Periodo; % Fatr in Hz
rATRI=abs(m); % ATRI in %