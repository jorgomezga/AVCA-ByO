function [rHomog] = homogen_ms (mInput,iLongMask)

% Gray-level homogeneity in the mInput image calculated using the Bhanu
% method, as stated in:
%
% [1] R. Peters and R. Strickland, “Image complexity metrics for automatic 
% target recognizers,” Autom. Target Recognizer Syst. Technol. Conf., 1990.
%
% INPUTS:
%       - mInput: image to detect contrast
%       - iLongMask: size of the mask (iLongMask x iLongMask)
%
% OUTPUT:
%       - rHomog: Global homogeneity value
%
% $Id: $

if nargin<2
    iLongMask=3;
end
if iLongMask<3, iLongMask=3; end

% average filter mask
mH=fspecial('average',iLongMask);

mU=(mInput-filter2(mH,mInput)).^2; %.^2 incluído en junio 2015
% Borders of the matrix are removed as homogeneity on borders can be
% confusing.
%iBorder=floor(iLongMask/2);
%mU=mU(1+iBorder:end-iBorder,1+iBorder:end-iBorder);
rHomog=mean(mean(mU));
        


