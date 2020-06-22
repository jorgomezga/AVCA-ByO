function A1 = LPCcoef1( vFrame )

% Calculates the first prediction coefficient considering a LPC analysis of order 12
% (coincides with the first sample of the cepstrum)
%
% Inputs:
%   vFrame               = Input vFrame
% Outputs:
%   vFo                  = Pitch value of the current frame

NFFT = 2^ceil( log2( length( vFrame ) ) );
vCepstrum = real( ifft( log( abs( fft(vFrame, NFFT))), NFFT) );

A1 = vCepstrum(1);


