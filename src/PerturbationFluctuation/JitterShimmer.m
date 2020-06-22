 function [vJitter, vShimmer] = JitterShimmer( vSignal, iFs )

% Calculation of the perturbation and fluctuation parameters
%
% Inputs:
%   vSignal              = Input signal
%   iFs                  = Sampling frequency
%
% Outputs:
%   vJitter              = Family of parameters based on jitter. Returns a
%                          vector with the following parameters:
%                           rJitta:   absolute value in microseconds
%                           rJitt:   relative value in % with respect to the average value of the period
%                           rRrRAP = rap( vPPS );
%                           rPPQ: pitch perturbation quotient in % 
%                           rSPPQ: smoothed pitch perturbation quotient in % (sPPQ).
%   vShimmer             = Family of parameters based on shimmer. Returns a
%                          vector with the following parameters:
%                           rShdB:     absolute shimmer in dB
%                           rShim:     relative shimmer in %
%                           rAPQ:      amplitude perturbation quotient % (APQ).
%                           rSAPQ:     smoothed amplitude perturbation quotient in % (sAPQ).

iLongitud  = length( vSignal );

% For the calculation of the disturbances and tremor some information is needed that
% is obtained from the pitch analysis
[vPPS, vPAS] = ParamPitch( vSignal, iFs, 1, iLongitud, 0 );

%% Pitch perturbation measures: jita, jitt, RAP, PPQ, sPPQ.
% rJitta:   absolute value in microseconds
% rJitt:   relative value in % with respect to the average value of the period
[rJitta, rJitt] = jitter( vPPS );
% rRAP: relative average disturbance of the pitch period in %
rRAP = rap( vPPS );
% rPPQ: pitch perturbation quotient in % 
rPPQ = ppq( vPPS );
% rSPPQ: smoothed pitch perturbation quotient in % (sPPQ).
rSPPQ=sppq( vPPS, 55 );
    
vJitter = [ rJitta, rJitt, rRAP, rPPQ, rSPPQ ];

%% Amplitude perturbation measures: shima, shim, APQ, sAPQ.
% rShdB:     absolute shimmer in dB
% rShim:     relative shimmer in %
[rShdB, rShim]=shimmer( vPAS );
% rAPQ:     amplitude perturbation quotient % (APQ).
rAPQ = apq( vPAS );
% rSAPQ:    amplitude perturbation quotient in % (sAPQ).
rSAPQ = sapq( vPAS, 55 );

vShimmer = [ rShdB, rShim, rAPQ, rSAPQ ];