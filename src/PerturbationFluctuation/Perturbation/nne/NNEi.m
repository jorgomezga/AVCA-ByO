function vNNEi=NNEi(s, Fs, fo, mBandas, bGraficar )

% Calculates the normalized noise energy (NNE) of a voice segment according to
% the Kasuya method [1]. Different bands of frequency can be used to perform the calculation.
%
% [1] H. Kasuya, S. Ogawa, K. Mashima and S. Ebihara. "Normalized noise energy
% as an acoustic measure to evaluate pathologic voice ", J.Acoust.Soc.Am.,
% vol. 80, no. 5, pp. 1329-1334, 1986.
%
% Input parameters:
%   s:    Row vector with the input voice segment. It has a duration long
%         enough to met the resolution in frequency required. It must be non-silent.
%   Fs:   the sampling frequency
%   fo:   the previously calculated pitch for the segment in Hz
%   m:    Matrix with the the lower and upper frequencies bands (in Hz)
%         which are used to calculate NNE is calculated. It should have two
%         rows and as many columns as bands
%   bGraficar: Plot the speech vs noise spectrum
% Output parameters:
% vNNEi: the normalized noise energy of the segment in the bands of
%           frequency requested.

if nargin < 5, bGraficar = 0;                           end
if nargin < 4, mBandas=[1000;fs/2];                     end
if nargin < 3, error( 'Not enough input parameters!' ); end

iNumBandas=size(mBandas,2);
N=length(s);

% Number of FTT coefficients
iNFFT=2^ceil(log2(N));

[S, S_peak, S_dip]=Noisei(s, Fs, fo);

%=======================================================================
%             NNE calculation for each frequency band
%=======================================================================
vNNEi = zeros(iNumBandas,1);
for b=1:iNumBandas
    fL = mBandas(1,b);     % Lower frequency (Hz)
    fH = mBandas(2,b);     % Higher frequency (Hz)
    
    nL = ceil(iNFFT*fL/Fs);   % Lower frequency in samples
    nH = floor(iNFFT*fH/Fs);  % Higher frequency in samples
    
    if nL == 0
        nL = 1;
    end
    
    E_dip = sum( S_dip(nL:nH) );
    E_peak = sum( S_peak(nL:nH) );
    E_ruido = E_dip + E_peak;
    E_tot = sum( S(nL:nH) );
    
    % The energy of the noise cannot be greater than the total energy of the signal.
    if E_tot < E_ruido
        E_ruido = E_tot;
    end
    
    vNNEi(b)=10*log10( eps+ ((E_ruido)/E_tot) );
    
    if bGraficar
        eje_f=((Fs/1000)/iNFFT).*(1:iNFFT/2);
        plot( eje_f(nL:nH), log( S(nL:nH) ), 'r' ),
        hold on,
        plot( eje_f(nL:nH), log( S_dip(nL:nH)+S_peak(nL:nH) ), 'g' )
        
        legend( {'Speech Spectrum', 'Noise Spectrum'} )
        title('log speech spectrum vs Noise (dB)');
    end
end