function vParameters = AdditiveNoise( vSignal, iFs, iFrameSize, iOverlap )

% Calculation of the additive noise parameters
%
% Inputs:
%   vFrame               = Input vFrame
%   iFs                  = Sampling frequency
%   iFrameSize           = Number of samples used to define the frame of
%                          analysis
%   iOverlap             = Number of samples for ovelapping windows
%
% Outputs:
%   vParameters          = Vector containing the following additive noise parameters 
%                           - vNNE: Normalized noise energy
%                           - vHNR: Harmonics to noise ratio - Yumoto
%                           - vGNE: Glottal-to-noise excitation ratio
%                           - vCHNR: Cespstrum harmonics-to-noise ratio

iLengthSignal = length(vSignal);
iNumVentanas = fix((iLengthSignal - iFrameSize + iOverlap)/(iOverlap));

% Normalized noise energy
vNNE  = NNE(vSignal, iFs, 1, iLengthSignal, iNumVentanas ); 
% Harmonics to noise ratio - Yumoto
vHNR = HNRYum(vSignal, iFs, 1, iLengthSignal, iNumVentanas); 
% Cespstrum harmonics-to-noise ratio
vCHNR = CHNR( vSignal, iFs, 1, iLengthSignal, iNumVentanas );

% Glottal-to-noise excitation ratio
vGNE = GNE(vSignal, iFs, 1, length(vSignal), iNumVentanas, 40e-3, 500, 100 ); 

vParameters = [vNNE, vHNR, vCHNR, vGNE];