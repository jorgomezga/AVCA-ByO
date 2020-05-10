function hmEnvolventes=CalculaEnvolventesHilbert( mTramas, iFs, iBW, iFShift, vBanda )

% Calculates the Hilbert envelopes in different frequency bands of each of the input voice
% frames. Frequency bands are determined by input parameters iFs, iBW, iFShift and vBanda.
% To calculate the envelope, the spectrum of the window is obtained,
% the spectrum divided into different frequency bands, the DFT is calculated
% as the inverse of each band and its module is taken.
%
% Input parameters:
%   mTramas:    An array that contains in its rows the temporary frames on
%               which Hilbert envelopes will be calculated.
%   iFs:        Sample rate (Hz) of the frames of mTramas
%   iBW:        Hilbert envelope bandwidth (Hz)
%   iFShift:    Shift or jump between frequency bands (Hz)
%   vBand:      Two-position vector with frequency range
%               (lower and upper, in Hz) for the calculation of the GNE
%
% Output parameters:
% 	hmEnvolventes:  Hyper-matrix of size (iNumBands x 2 * iWidthBand x
%                   iNumVentanas) with the calculated Hilbert envelopes.
%                   The 1st index indicates the frequency band and the 3rd
%                   the time window to which it belongs.

if nargin < 5, vBanda = [0 5000]; end
if nargin < 4, iFShift=100; end
if nargin < 3, iBW=3000; end
if nargin < 2, iFs=10000; end

% Number of temporary windows and samples of each window.
[iNumVentanas, iNumMuestras] = size( mTramas );

% Number of FFT coefficients
iNFFT = 2^nextpow2( iNumMuestras );

% Number of discrete frequency samples corresponding to iBW.
iAnchoBanda = ceil( iBW*(iNFFT/iFs) ); 

% Number of discrete frequency samples corresponding to iFShift.
iDesplazamiento = ceil( iFShift*(iNFFT/iFs) );

% GNE calculation start frequency (in samples).
iInicio = ceil( vBanda(1)*(iNFFT/iFs) );

% GNE calculation end frequency (in samples).
iFinal = min(iNFFT/2, ceil( vBanda(2)*(iNFFT/iFs) ));

% Number of frequency bands
iNumBandas = floor( ((iFinal-iInicio)-iAnchoBanda)/(iDesplazamiento) ) + 1;

% Create frequnecy bands
vVentana = hann( iAnchoBanda )';
mVentanas = zeros( iNumBandas, iNFFT );
for b=1:iNumBandas
    mVentanas( b, 1+iInicio+(iDesplazamiento*(b-1)):iInicio+(iDesplazamiento*(b-1))+iAnchoBanda ) = vVentana;
end

% FFT of temporal windows
mSenalesFFT = fft( mTramas', iNFFT ).';

% Inicialization of the hypermatrix
hmEnvolventes = zeros(iNumBandas, iNumMuestras, iNumVentanas);

for v=1:iNumVentanas
    for b=1:iNumBandas
        % Frequency windowing to obtain a band pass signal.
        vTemp = mSenalesFFT(v,:).*mVentanas(b,:);
        % We obtain the b-th Hilbert envelope.
        vTemp = abs( ifft( vTemp, iNFFT ) );
        hmEnvolventes(b,:,v) = vTemp( 1:iNumMuestras );
    end
end