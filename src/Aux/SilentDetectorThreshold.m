function bSound = SilentDetectorThreshold( vFrame )

% Calculates if the vFrame is silent or not according to the energy of
% the segment using thresholds
%
% Inputs:
%   vFrame               = Input vFrame
% Outputs:
%   bSound               = Boolean indicating if the frame is silent or not

% Zero-crossing rate of the segment
rZCR = ZeroCrossingRate( vFrame );
% Log-Energy of the segment
logEnergy = LogEnergy( vFrame );

% Thresholds for voice energy
THRESH_VOICE_ENER=-40;
THRESH_VOICE_ZCR=0.3;

bSound = false;

% Si el segmento supera un umbral de energ�a o de TCC se considera que es voz y se
% analiza su sonoridad, sino se termina y se devuelve S=0
if logEnergy > THRESH_VOICE_ENER...
        || rZCR > THRESH_VOICE_ZCR
    % If the first condition is met then, other conditions are considered to 
    % analyze if the segment is voiced or not:
    
    % Thresholds for the voiced/unvoiced decision based on energy thresholds
    THRESH_LOGENER=-25;    
    % zero-crossing rate threshold
    THRESH_ZCR=0.07;
    
    % First coefficient of the autocorrelation function
    rC1 = AC1( vFrame );
    % Threshold
    THRESH_C1=0.8;
    
    % First LPC coefficient
    rA1 = LPCcoef1( vFrame );
    % Threshold
    THRESH_A1=-1;
    
    % Energy of the normalized prediction error
    iPredictionOrder = 12;
    rEp = ErrorPredNorm( vFrame, iPredictionOrder );
    % Threshold
    THRESH_PRED_ERR=10;
   
    % Count the number of conditions met
    iNumConditions = (logEnergy > THRESH_LOGENER) + (rZCR < THRESH_ZCR) +...
            (rC1 > THRESH_C1) + (rA1 < THRESH_A1) + (rEp > THRESH_PRED_ERR);
    
    if iNumConditions>=2
        bSound=true;
    end
else  
    % If the frame does not fullfil the THRESH_VOICE_ENER and THRESH_VOICE_ZCR
    % threshold, then the segment is silent
    bSound = false;    
end