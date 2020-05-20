function vmOut=sintonizaEM(vmMS,vMf,vAf,iFmodMax,iFaMin,iFaMax)

% Limits Modulation Spectra (vmMS) to the specified modulation (iFmodMax)
% and acoustic (iFaMin, iFaMax) ranges.
%
% INPUTS:
%           vmMS - matrix vector containing the modulation spectrum.
%           vMf - vector containing the modulation frequency bins
%           vAf - vector containing the acoustic frequency bins
%           iFmodMax - maximum modulation frequency
%           iFaMin - minimum acoustic frequency
%           iFaMax - maximum acoustic frequency
%
% OUTPUTS:
%           vmOut - matrix vector containing the output modulation
%           spectrum, with new acoustic and modulation frequency limits
%

% Modulation frequency filtering
if iFmodMax<vMf(end)
    sSpecopt=max(find(vMf<=iFmodMax));
    vmMS=vmMS(:,1:sSpecopt,:);
end

% Acoustic frequency filtering

[~,iAcmin]=min((abs(vAf-iFaMin)));
[~,iAcmax]=min((abs(vAf-iFaMax)));

if iAcmin==1, iAcmin=0;end

vmOut=vmMS(iAcmin+1:iAcmax,:,:);

