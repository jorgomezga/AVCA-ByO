function [ mFeatures, caFeatureNames, eMS] = MS_features( vSignal, iFs, sTipo, iSubMod, iFrame, iShift, eOptions, iVerbosity)

% Generates the Modulation Spectra (MS) of an input signal and its
% MS features as detailed in refs [1] and [2]
% INPUTS:
%       vSignal: vector containing the input signal (mono)
%       iFs - The sampling rate of vSignal, in Hz.
%       sTipo    any sensible combination of the following:
%               'R'  rectangular window in time domain
%               'N'  Hanning window in time domain
%               'M'  Hamming window in time domain (default)
%               't'  triangular shaped filters in mel domain (default)
%               'n'  hanning shaped filters in mel domain
%               'm'  hamming shaped filters in mel domain
%               'C' calculate Modulation Spectra Centroids
%               'Y' calculate Modulation Spectra Dynamic Range for every
%                   Modulation Frequency Band
%               'f' normalize Centroids respect to pitch
%               'F' normalize Centroids respect to the 0Hz centroid
%               'O' calculate "Low Modulation Ratio" in Modulation Spectra.
%               'P' calculate Contrast, and homogeneity in Modulation Spectra modulus
%                   and phase
%               'A' Calculate Modulation Spectra Histogram features.
%        iSubMod - Number of centroids
%        iFrame - framelength in samples (minimum has to be equivalent to
%        40 ms. Otherwise MS will not contain relevant information about
%        low modulation frequency.
%        iShift - frame shift length in samples
%        eOptions - Structure containing other auxiliar and secondary
%        inputs, such as:
%               - iFmodMax - maximum modulation frequency
%               - vFlmr - max frequency used in LMR calculation
%               - EMiSubbands - desired number of acoustic subbands
%               - EMsSpecopt - desired number of modulation bands
%               - EMvFlmr - High bound of the band used to calculate LMR
%               (25 Hz by default)
%               - EMsDemod - A data structure containing demodulation options. This can
%               be a string indicating the demodulation method, or
%               alternatively a cell array specifying parameter values in
%               the fashion of the MODDECOMP... functions. The default
%               setting is {'cog', 0.1, 0.05}.
%                   {'HILB'}
%                   {'COG', <carrwin>, <carrwinhop>}
%                       carrwin - seconds.
%                       carrwinhop - seconds.
%                   {'HARM', <numharmonics>, <voicingsens>, <F0smoothness>}
%                       numharmonics - a positive integer.
%                        voicingsens - a decimal value (0 to 1).
%                       F0smoothness - a positive integer.
%                   {'HARMCOG', <carrwin>, <carrwinhop>, ...
%                               <numharmonics>, <voicingsens>, <F0smoothness>}
%   	 iVerbosity - if iVerbosity==1, the function displays texts related
%   	 to the calculation process in screen
%
% OUTPUTS:
%       mFeatures - Matrix containing the MS features. Each row will
%       be referred to a frame and, depending on the values of sTipo,
%       will contain the features:
%       [Centroids, Dynamic range, LMR, Modulus contrast, Modulus homogeneity,
%           Phase contrast, Phase homogeneity, CIL, PALA, RALA, RALP25,
%           RALP75, Ralp95]
%       eMS - Structure containing:
%           vmMS - matrix vector containing the MS of each signal frame.
%           vMf - modulation frequencies vector
%           vAf - acoustic frequencies vector
%           data - data structure from modspectrum function
%       caFeatureNames - Includes the names of the calculated features. Has
%       the same length as the number of columns of mFeatures.
%
% REFERENCES:
% [1] Moro-Velázquez, L., Gómez-García, J. A., Godino-Llorente, J. I., &
% Andrade-Miranda, G. (2015). Modulation spectra morphological parameters:
% a new method to assess voice pathologies according to the grbas scale.
% BioMed research international, 2015.
%
% [2] Moro-Velázquez, L., Gómez-García, J. A., & Godino-Llorente, J. I.
% (2016). Voice pathology detection using modulation spectrum-optimized
% metrics. Frontiers in bioengineering and biotechnology, 4, 1.

if nargin<2, error('Please, include a sampling frequency' ); end
if nargin<3, sTipo='MCYFOPA'; end % All features are calculated, using Hamming windowing
if nargin<4, iSubMod=26; end
if nargin<5, iFrame=pow2(floor(log2(0.150*iFs))); end % Default: 150 ms
if nargin<6, iShift=floor(iFrame/2); end
if nargin<7, eOptions=[]; end
if nargin<8, iVerbosity=0; end

if isfield ( eOptions, 'EMiSubbands' ), iSubbands=eOptions.EMiSubbands;
else, iSubbands= 128; end

if isfield ( eOptions, 'EMsSpecopt' ), sSpecopt=eOptions.EMsSpecopt;
else, sSpecopt= 1024; end

if isfield (eOptions, 'EMsDemod' ), sDemod=eOptions.EMsDemod;
else, sDemod= 'HILB'; end

%iFil: input signal length
[iFil]=size(vSignal,1);

if iFrame>iFil
    iFrame=iFil;
end

if iFrame<round(iFs*0.04)
    iFrame=round(iFs*0.04);
    warning(['iFrame must be longer than 40 ms of samples. Changing iFrame to: ' num2str(iFrame)])
end

%% Framing plus windowing
[mSignal]=enframing_avca( vSignal, iFrame, iShift, sTipo );

iNumtramas=size(mSignal,1); % Number of frames

if length(iSubbands)>1
    vBandas=iSubbands; % vBandas contains subband acoustic frequency boundaries
    iNumBands=length(vBandas)-1;
else
    vBandas=(iFs/2)/(iSubbands-1); % vBandas contains acoustic bandwidth for uniform-width
    %               subbands
    
    iNumBands=iSubbands;
end

% Output matrix vector
vmMS=zeros(iNumBands, sSpecopt,iNumtramas);
eOptions.iShift=iShift;

for j=1:iNumtramas
    
    [vmMS(:,:,j),vMf,vAf,data]=...
        modspectrum(mSignal(j,:),iFs, sDemod, vBandas, sSpecopt);
    % vMf y vAf are respectively the same for every frame.
end

eMS=struct('vmMS', vmMS, 'vMf', vMf, 'vAf', vAf, 'data', data);

if iVerbosity
    disp('Modulation spectrum feature calculation')
end

% Once vmMS is calculated, MS features can be calculated:
[ mFeatures, caFeatureNames ] = EM_ePar( eMS, vSignal, sTipo, iSubMod, eOptions, iVerbosity);

if nargout<3
    eMS=[];
end
