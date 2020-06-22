function HNRi=CHNRKromi( vFrame, iFs )

% Calculates the harmonic-to-noise ratio (HNR) of a voice segment, using
% the method used by Guus de Krom.
%
% Input parameters
%       vFrame: the voice frame. It is supposed to be at least as long
%           of the longest period in voice (15ms). Another option is to
%           use the previously calculated pitch; since it is enough
%           to have more than one signal period, having to search for
%           the main rahmonic in the neighourhood of the pitch period.
%           NOTE: HNR depends on the size of the analysis window.
%       iFs: the sampling frequency
%
% Output parameters
%       HNRi:     the value of the HN

if nargin < 2, error( 'Not enough input parameters!' ); end

% Check that the vector is of type row 
if ~isvector( vFrame )
    error( 'El parámetro vSignal no es un vector columna' ); 
elseif size( vFrame, 2 ) == 1
    vFrame=vFrame'; 
end 

% find and fsup are the lower and upper frequencies (in Hz) of the
% frequencies range at which the HNR is calculated.
% Frequency range by default in Krom algorithm: 60-2000 Hz
finf = 60; 
fsup = 2000;

N = length( vFrame );

% Number of FFT coefficients
NFFT=2^( ceil( log2(N)));
% If N< NFFT -> Zero padding
if N<NFFT
    vFrame=[vFrame zeros(1,NFFT-N)];
end

% Hanning windowing
sw=vFrame.*hanning(NFFT)';

% Log Spectrum: O(f).
S    = fft(sw,NFFT);
magS = abs(S);
log_espec=log(magS);

% Cepstrum
cepstrum=real(ifft(log_espec));

% Only the real part is kept -> Half of the signal is discarded
magS=magS(1:NFFT/2);
log_espec=log_espec(1:NFFT/2);
%cepstrum=cepstrum(1:NFFT/2);

% Location of the main rahmonic, corresponding to the pitch;
% this value is looked for between 500 Hz (2 ms) and 66.6 Hz (15 ms).
mueMin= round(0.002*iFs);
mueMax= min( round(0.015*iFs), NFFT/2 );

[~,pos] = max(cepstrum(mueMin:mueMax));
rahmain = pos+(mueMin-1);

% Location of the next rahmonics; looking in multiples of
% pitch, within a 20% margin. Their positions are stored in the variable
% rahmonicos. The search procedure is carried out only in the first half 
% of the cepstrum as it is real.
margen=round(0.2*rahmain);
rahmonicos(1)=rahmain;
i=2;
while i*rahmain < NFFT/2
    mueMin=i*rahmain - margen;
    mueMax=min(i*rahmain + margen,NFFT/2);
    [~,pos]=max(cepstrum(mueMin:mueMax));
    rahmonicos(i)=pos+(mueMin-1);
    i=i+1;
end

% Bandwith of the rahmonics
iNumRahmonicos = length(rahmonicos);
B = zeros(1, iNumRahmonicos);
for i=1:iNumRahmonicos
    d=rahmonicos(i);
    while (d<NFFT/2) && (cepstrum(d)>0)
        d=d+1;
    end
    B(i)=2*(d-rahmonicos(i));
end

% Calculation of the cepstrum without the part corresponding to the rahmonics
% (filtered in the cepstral domain: "comb-liftering").
no_rahm=cepstrum;
for i=1:iNumRahmonicos
    for d=rahmonicos(i)-B(i)/2 : rahmonicos(i)+B(i)/2
        if d>0
            no_rahm(d)=0;
        end
    end
end

% The last half of the "liftered" cepstrum is obtained from the first
% half of the frame since it is symmetrical. With a symmetric cepstrum c(r) and r=0...NFFT-1,
% we have that c (r) = c (NFFT-r). In Matlab: c(r+1)=c(NFFT-r + 1).
for r=1:NFFT/2-1
    % Nota: no_rahm(NFFT/2+1) no lleva pareja.
    no_rahm(NFFT-r+1)=no_rahm(r+1);
end

% Calculation of the noise spectrum, Nap.
Nap=real(fft(no_rahm));

% Only half of the samples are representative, as it is real no_rahm.
Nap=Nap(1:NFFT/2);

% Calculation of the approximation to the harmonic spectrum, Haap.
Haap=log_espec-Nap;

% Calculation of the correction terms Bd(f), for the baseline correction
% which removes the spectral envelope information present in Nap(f),
% subsequently obtaining the noise spectrum N (f).

% ko is the frequency sample corresponding to pitch fo.
ko = NFFT/rahmain;

% The maximums in Haap corresponding to the harmonics will have values 
% close to the multiples of ko within a 30% margin of said multiples
% The result is stored in the variable k(i).
margen = round(0.3*ko);
k = zeros(1,length(Haap)); 
i=1;
while round(i*ko) < length(Haap)
    mueMin = max( round(i*ko) - margen, 1 );
    mueMax = min( round(i*ko) + margen, length(Haap) );
    [~,pos] = max(Haap(mueMin:mueMax));
    k(i) = pos+(mueMin-1);
    i = i+1;
end
k(i:end) = []; % Clipped to the real length

% The minimum (in absolute value) of Haap are calculated between successive
% harmonics, obtaining a series of corrective terms Bd (f) such
% that: Bd (f) = | min [Haap (f)] |, where f varies as (n-1) fo <f <= n fo
Nk = length(k);
M = zeros(1, Nk+1);
Bd = zeros(1, length(Haap));

% The first absolute minimum is calculated between the first Haap sample and
% its first harmonic.
M(1) = abs(min(Haap(1:k(1))));
Bd(1:k(1))=M(1);

% The absolute minimums between successive harmonics are obtained.
for i=2:Nk
    M(i)=abs(min(Haap( k(i-1):k(i))));
    Bd(k(i-1):k(i))=M(i);
end

% The last absolute minimum is calculated between the final and the last harmonic.
% sample of Haap
M(Nk+1)=abs(min(Haap(k(Nk):length(Haap))));
Bd(k(Nk):length(Haap))=M(Nk+1);

% Noise spectrum:
NoiseSpectrum=Nap-Bd;
magN=exp(NoiseSpectrum);%10.^N;

% kind is the sample corresponding to the lower frequency of the
% margin. If we analyze from the continuous component (f = 0) we take kinf = 1.
kinf=round(finf*(NFFT/iFs));
if kinf==0
    kinf=1;
end

% kups is the number of samples corresponding to the upper frequency of the
% margin. The rounding out of the admissible limit is prevented: NFFT / 2
ksup=round(fsup*(NFFT/iFs));
if ksup>NFFT/2
    ksup=NFFT/2;
end

% HNR calculation (according to the Guus algorithm) at the provided
% frequency range
signalPot=0;
noisePot=0;
for i=kinf:ksup
    signalPot=signalPot+(magS(i)^2);
    noisePot=noisePot+(magN(i)^2);
end
if noisePot>0
    HNRi=10*log10(signalPot/noisePot);
else
    HNRi=0;
end

% Plots
if nargout == 0
    eje_f=((iFs/1000)/NFFT).*(1:NFFT/2);
    eje_ceps=1:NFFT;
    figure(1); plot(eje_f,log10(exp(log_espec)),'r',eje_f,log10(exp(Nap)),'g');
    title('O(f) vs Nap(f)');
    figure(2); plot(eje_f,log10(exp(log_espec)),'b',eje_f,log10(exp(Haap)),'g');
    title('log espectro vs Haap en dB');
    figure(3); plot(eje_f,log10(exp(log_espec)),'r',eje_f,log10(exp(Bd)),'g');
    title('log espectro vs Bd en dB');
    figure(4); plot(eje_f,log10(exp(log_espec)),'r',eje_f,log10(exp(NoiseSpectrum)),'g');
    legend( {'Speech Spectrum', 'Noise Spectrum'} )
    title('log espectro vs Noise en dB');
    figure(5); plot(eje_ceps,cepstrum,'r',eje_ceps,no_rahm,'g');
    title('liftrado del cepstrum');
    figure(6); plot(eje_f,log10(exp(Haap)),'r',eje_f,log10(exp(Bd)),'g');
    title('Haap vs Bd en dB');
end
