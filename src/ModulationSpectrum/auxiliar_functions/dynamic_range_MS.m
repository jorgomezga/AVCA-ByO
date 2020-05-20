function [vDR]=dynamic_range_MS(mModSpectra)

% function [vDR]=dynamic_range_MS(vModSpectra)
%
% Calculates Dynamic Range for every modulation frequency band included in
% input matrix mModSpectra.
%
% INPUTS:
%       - mModSpectra: matrix including modulation spectra energy values
%       (dB).
%
% OUTPUTS:
%       - vDR: vector including dynamic range of input matrix (dB)
%

vDR=max(mModSpectra,[],1)-min(mModSpectra,[],1);



