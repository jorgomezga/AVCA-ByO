function [vAmax, vPosAmax, mPicosPos]=PicosMayores( vSignal, iFs, iInicio, iFinal, lMetodo )

% Calculate the highest energy peaks, cycle by cycle of the input voice
%
% Input parameters
%   vSignal:       is a column vector containing the complete speech signal
%   iFs:           is the sampling frequency
%   iInicio:       is the first sample of the section to be analyzed
%   iFinal:        is the last sample of the section to be analyzed
%   lMethod:       indicates the method to be used:
%                   lMethod = 0: Boyanov
%                   lMethod = 1: Kasuya
%
% Output parameters
%   vAmax:         amplitudes of the highest energy peaks, cycle by cycle, of a speech frame;
%   vPosAmax:      positions of the highest energy peaks, cycle by cycle, of a voice frame;
%   mPicosPos:     Peak signal (of the same duration as vSignal) with the peaks found.

% Check that the vector is of type column
if ~isvector( vSignal )
    error( 'vSignal is not a vector!' );
elseif size( vSignal, 2 ) ~= 1
    vSignal=vSignal';
end

if nargin < 5, lMetodo=0;                   end
if nargin < 4, iFinal=length( vSignal );    end
if nargin < 3, iInicio=1;                   end
if nargin < 2, error( 'Not enough input parameters!' ); end

switch lMetodo
    case 0 % (Boyanov)
        
%         % iFrame = 40ms
%         iFrame=fix(0.034*iFs);
%         % iDesplazamiento = 20ms (overlpa of 50%).
%         iDesplazamiento=iFrame;%/2;
%         % Total number of windows
%         iNumVent=ceil( length( vSignal )/iDesplazamiento );
%         [~, vF0] = CalculatePitch( vSignal, iFs, iNumVent, 'Boyanov' );
%         vT = round( iFs./vF0 );

        vT = PitchBoy( vSignal, iFs, iInicio, iFinal );
        
        % Postive peaks
        [vApos, vPpos, vPicosPos]=Picos( vSignal, iFs, iInicio, iFinal, vT, 1);
        
        % Negative peaks
        [vAneg, vNneg, vPicosNeg]=Picos( vSignal, iFs, iInicio, iFinal, vT, -1);
        
    case 1 % (kasuya)
        
%         % iFrame = 40ms
%         iFrame=fix(0.04*iFs);
%         % iDesplazamiento = 20ms (overlap of 50%).
%         iDesplazamiento=iFrame/2;
%         % Total number of windows
%         iNumVent=ceil( length( vSignal )/iDesplazamiento );
%         [~, vF0] = CalculatePitch( vSignal, iFs, iNumVent, 'Kasuya' );
%         vT = round( iFs./vF0 );
        
        % The sequence of pitch periods is obtained by the Kasuya (Feijoo) method, 
        % with 40ms windows and 50% overlap.
        vT=PitchKas( vSignal, iFs, iInicio, iFinal );        
        
        % Null values that exist at the beginning and ending of the frame are removed if the
        % first sample of signal s or last sample are included in the analysis.
        % (The value 0 is arbitrarily assigned if the analysis window is left out
        % the voice signal. Remember that windows are taken centered at the point of
        % analysis).
        vT=QuitaExtremos( vT );
        if isempty(vT)
            vAmax= [];
            vPosAmax=[];
            mPicosPos=[];
            warning( 'It is not possible to calculate the pitch period of the input signal, returning 0 value' )
            return
        end
        
        % The positive peaks of each period of the considered segment are obtained.
        [vApos, vPpos, vPicosPos]=PicosKas( vSignal, iFs, iInicio, iFinal, vT, 1);
        
        % The negative peaks of each period of the considered segment are obtained.
        [vAneg, vNneg, vPicosNeg]=PicosKas( vSignal, iFs, iInicio, iFinal, vT, -1);
end

% We will keep the peaks that have the highest mean absolute value (highest energy)
if ~isempty( vApos )
    rEpos=sum( vApos )/length( vApos );
else
    rEpos=0;
end

if ~isempty( vAneg>0 )
    rEneg=abs(sum( vAneg ))/length( vAneg );
else
    rEneg=0;
end

if rEpos>= rEneg
    vAmax=vApos;
    vPosAmax=vPpos;
    mPicosPos=vPicosPos;
else
    vAmax=abs( vAneg );
    vPosAmax=vNneg;
    mPicosPos=vPicosNeg;
end

if nargout ==0
    vEje=( iInicio:iFinal )/iFs;
    plot( vEje, vSignal( iInicio:iFinal ),'r', vEje, mPicosPos( iInicio:iFinal ),'g');
    title('Major peaks');
end