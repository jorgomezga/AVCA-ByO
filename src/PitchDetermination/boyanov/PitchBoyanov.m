function iFo = PitchBoyanov( vFrame, iFs, iT0 )

% Calculates the pitch of a given frame using the method described in:
% B. Boyanov, S. Hadjitodorov, B. Teston, D. Doskov. Robust hybrid pitch
% detector for pathologic voice analysis. Larynx 97, Jun 1997, Marseille, France. pp.55-58.
% According to the temporal and cesptral method
%
% Inputs:
%   vFrame               = Input vFrame
%   iFs                  = Sampling frequency
%   iT0                  = Previous calculation of the pitch period
% Outputs:
%   vFo                  = Pitch value of the current frame

if nargin < 2, error( 'Not enough input parameters!' ); end
if nargin < 3, iT0 = 0;             end

% Silence/Non-silence decision making based on thresholding
bSound = SilentDetectorThreshold( vFrame );

% Additional voiced/unvoiced threshold
VOICED_THRESHOLD=0.3;

if bSound
    % Central clipping with dynamic adaptation of clipping thresholds
    % eliminate noise and formant components difficulting F0 calculation.
    vProc = DynamicClip( vFrame );
else
%     iFo = 0;
    iFo = NaN;
    return
end

% The autocorrelation function of the clipped signal (temporal) or the
% liftered signal (cepstral) is computed
vAutocorr = xcorr( vProc );

% The pitch is obtained from the autocorrelation function
iFo = PitchDetecBoy( vAutocorr, iFs, VOICED_THRESHOLD, iT0 );