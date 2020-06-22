function [S, S_peak, S_dip, len_rdos]=Noisei(s,Fs,fo)

% This function allows calculating the energy of the noise present in the segment,
% when separating the peaks of the spectrum of the valleys.
%
% Input parameters
%   s:  row vector containing the speech segment (NNE) or its Hilbert transformation
%       It has to be long enough to have the necessary frequency resolution. 
%       It must be a sound segment and the average value must have been subtracted
%       and Hamming windowed.
%   Fs: sampling frequency
%   fo: the previously calculated pitch for the segment in Hz
%
% Output parameters
%   S:      the total power spectrum of the voice segment
% S_peak:   the noise power spectrum in the peak regions
% S_dip:    the noise power spectrum in the valley regions
% len_rdos: length of the power spectra

if nargin < 3
    error('Not enough input parameters'); 
end

% Si el vector no es de tipo columna se supone que se ha pasado como fila 
if ~isvector( s )
    error( 'Parameter s is not a vector' ); 
elseif size( s, 2 ) ~= 1 
    s=s'; 
end

M=length(s);

% This value is related to the width of the Hamming window that is
% used to get the width of the harmonics
Bw=2;

% The number of samples for the DFT is the value, power of two, greater than the
% segment length
NFFT=2^ceil(log2(M));

% Half the hamming window
semi_ham=round(Bw*NFFT/M);    

%=======================================================================
%						Calculation of the power spectrum
%=======================================================================
S=abs(fft(s,NFFT)).^2;
% If the input signal is a voice (real), the fft is symmetric and is sufficient
% with the just the first half; if the input signal is the Hilbert transformation
% of the voice signal (complex), the first half of the fft has twice the
% energy, and the second half zero -> In any case we are only interested in the first half.
% NOTE: the spectral resolution with the Hilbert transformation is the same
S=S(1:NFFT/2+1);
len_rdos=length(S);

%=======================================================================
%		Calculation of the power spectrum peaks (harmonics)
%=======================================================================
ko=round(NFFT*fo/Fs);
margen=round(0.4*ko);
i=1;
while i*ko+semi_ham < length(S)
   mueMin=i*ko - margen;
   mueMax=min( i*ko + margen, length(S) );
   [m,pos]=max(S(mueMin:mueMax));
   k(i)=pos+(mueMin-1);
   i=i+1;
end
N_picos=length(k);

%=======================================================================
%				Calculation of the ends of the peak regions
%=======================================================================
inf_ham=k-semi_ham;  
% Check that inf_ham is larger than 0.
inf_ham( inf_ham<=0 ) = 1; 

sup_ham=k+semi_ham;
% We will include the NFFT / 2 + 1 value in inf_ham, for the last DIP region.
inf_ham=[inf_ham,NFFT/2+1];
% We include in sup_ham the value semi_ham + 1 that is necessary for the first
% DIP region.
sup_ham=[semi_ham+1,sup_ham];  

%=======================================================================
% Calculation of the noise power spectrum in the DIP regions (valleys)
%=======================================================================
S_dip=S;

for i=1:N_picos
   S_dip(inf_ham(i):sup_ham(i+1))=0;
end

S_dip(1:semi_ham+1)=0; 

%=======================================================================
% Calculation of the noise power spectrum in the PICOS regions		
%=======================================================================
S_peak=zeros(size(S));
for i=1:N_picos
   
   if (sup_ham(i)+1<=inf_ham(i)-1) && (sup_ham(i+1)+1<=inf_ham(i+1)-1)
      medioL=mean(S(sup_ham(i)+1:inf_ham(i)-1));
      medioR=mean(S(sup_ham(i+1)+1:inf_ham(i+1)-1));
      Medio=0.5*(medioL+medioR);
   elseif (sup_ham(i)+1<=inf_ham(i)-1) && (sup_ham(i+1)+1>inf_ham(i+1)-1)
      medioL=mean(S(sup_ham(i)+1:inf_ham(i)-1));
      Medio=medioL;
   elseif (sup_ham(i)+1>inf_ham(i)-1) && (sup_ham(i+1)+1<=inf_ham(i+1)-1) 
      medioR=mean(S(sup_ham(i+1)+1:inf_ham(i+1)-1));
      Medio=medioR;
   else
      Medio=0;
   end
     
   S_peak(inf_ham(i):sup_ham(i+1))=Medio;
    
end