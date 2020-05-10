function [rElog, rE] = LogEnergy( vFrame )

% Calculates the energy of a frame in logarithmic and linear scale, 
% normalized by the number of samples
%
% Inputs:
%   vFrame     = Input vFrame
% Outputs:
%   rElog      = Energy in logarithmic scale
%   rE         = Energy in linear scale. 

rE    = (1/length( vFrame ))*sum( vFrame.^2 );
rElog = 10*log10( eps+rE );