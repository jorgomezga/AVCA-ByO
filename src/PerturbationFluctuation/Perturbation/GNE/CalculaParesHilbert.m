function mParejas = CalculaParesHilbert( iNumMuestras, iFs, iBW, iFshift, vBanda )

% Calculates the center frequencies for which the cross correlations of Hilbert envelopes are to be computed.
% Returns an array with two columns, on which one of the
% indices takes values between 1 and iNumBands, matching those frequency bands (windows) 
% whose center frequencies meet the next criteria: their difference is greater than or 
% equal to half the bandwidth or frequency window
%
% Input parameters:
%   iNumSamples:      Number of samples in each time window on which
%                     the Hilbert envelopes are to be computed.
%   iFs:              Sample rate (Hz) of the windows.
%   iBW:              Bandwidth of the Hilbert envelopes (Hz).
%   iFshift:          Shift or jump between frequency bands (Hz).
%   vBanda:            Two-position vector with frequency margins
%                     (lower and upper, in Hz) for the calculation of the GNE
%
% Output parameters:
%   mParejas:        Matrix with the correlation pairs and a two columns that store the
%                    indices (from 1 to iNumBands) relating the correlation pairs.
%
% See GNE.m

if nargin < 5, vBanda = [0 5000]; end
if nargin < 4, iFshift=100; end
if nargin < 3, iBW=3000; end
if nargin < 2, iFs=10000; end 

% Number of FFT coefficients
iNFFT = 2^nextpow2( iNumMuestras );

% Número de muestras en frecuencia discreta correspondientes a iBW.
iAnchoBanda = ceil( iBW*(iNFFT/iFs) ); 
                
% Number of discrete frequency samples corresponding to iFshift.
iDesplazamiento = ceil( iFshift*(iNFFT/iFs) );

% GNE calculation start frequency (in samples).
iInicio = ceil( vBanda(1)*(iNFFT/iFs) );

% GNE calculation end frequency (in samples).
iFinal = min(iNFFT/2, ceil( vBanda(2)*(iNFFT/iFs) ));

% Number of frequency bands
iNumBandas = floor( ((iFinal-iInicio)-iAnchoBanda)/(iDesplazamiento) ) + 1;

% Number of samples used in discrete frequency.
iLongEspectro = (iNumBandas-1)*iDesplazamiento+iAnchoBanda; 

% Vector with the iNum Bands starting points of the frequency bands.
vIndicesInicioBandas = (1:iDesplazamiento:iLongEspectro-iAnchoBanda+1);
vIndicesInicioBandas = vIndicesInicioBandas + iInicio -1;
                                                                                      
% Vector with the iNumBand indices of the center frequencies.
vIndicesFrecCentrales = vIndicesInicioBandas+ceil((iAnchoBanda+1)/2)-1;

% Vector with iNumBandas center frequencies in Hz.
vFrecCentrales = vIndicesFrecCentrales*(iFs/iNFFT); 

% Condition that two bands must meet to calculate their correlation.
rCondicion = iBW/2;

mParejas = zeros( iNumBandas*iNumBandas, 2 );
k = 1;
for i=1:iNumBandas
    for j=i:iNumBandas
        if abs( vFrecCentrales(i)-vFrecCentrales(j) ) >= rCondicion
            mParejas(k,1) = i;
            mParejas(k,2) = j;
            k = k+1;
        end
    end
end

mParejas(k:end,:) = [];