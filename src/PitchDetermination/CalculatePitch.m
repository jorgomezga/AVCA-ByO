function [iF0Out, vF0] = CalculatePitch( vSignal, iFs, iNumberOfFrames, sType )

% Calculates the pitch of a given frame using different pitch
% implementations
%
% Inputs:
%   vSignal              = Input signal
%   iFs                  = Sampling frequency
%   iNumberOfFrames      = desired number of frames. The overlapping is
%                          varied to ensure that the number of frames is
%                          achieved. The default value, iDesiredNumberOfFrames, 
%                          is twice the number of the frames that will result
%                          if the window size was fit to 40e-3. This emulates an
%                          overlap of 50%.
%
%   sType  =    'Boyanov': Dynamic clipping in the temporal domain
%               according to:
%                   B. Boyanov, S. Hadjitodorov, B. Teston, D. Doskov. Robust hybrid pitch
%                   detector for pathologic voice analysis. Larynx 97, Jun 1997, Marseille, France. pp.55-58.
%
%               'Rabiner_StaticClip': Static clipping -> Signal is clipped using
%                       static thresholds
%               'Rabiner_LPC'
%               'Rabiner_Cesptral'
%               'Kasuya' according to: 
%                   Kasuya, H., Ogawa, S., Kikuchi, Y., & Ebihara, S. (1986). 
%                   An acoustic analysis of pathological voice and its application 
%                   to the evaluation of laryngeal pathology. 
%                   Speech communication, 5(2), 171-181.
% Outputs:
%   vFo  = Pitch value

if nargin < 2, error( 'Not enough input parameters!' ); end
if nargin < 3, iNumberOfFrames = fix( 2*length( vSignal )/(iFs*40e-3) );    end
if nargin < 4, sType = 'Boyanov_Temporal';      end

vSignal = normalize( vSignal, 'zscore' );

% Length of the input signal
iLengthSignal = length( vSignal );

% The overlap is calculated in accordance to the desired number of frames
iOverlap = iLengthSignal/iNumberOfFrames;

switch (sType)
    case {'Boyanov'}
        % Boyanov -> Dynamic clipping: Varying iFrameSize
        iFrameSize = round(34e-3*iFs);
        iNumberTo=3;
    case 'Rabiner_StaticClip'
        % Rabiner -> Static clipping: Varying iFrameSize
        iFrameSize = round(34e-3*iFs);
        iNumberTo=5;
    case 'Rabiner_LPC'
        % 'Rabiner_LPC': Fixed iFrameSize
        iFrameSize = round(60e-3*iFs);
    case 'Rabiner_Cepstral'
        % Rabiner -> Cepstral: Fixed iFrameSize
        iFrameSize = round(40e-3*iFs);
    case 'Kasuya'
        % Kasuya -> :Fixed iFrameSize
        % For this particular case, the framesize and the overlap are already fixed
        iFrameSize = round(40e-3*iFs);
        iOverlap   = round(20e-3*iFs);
        
        nli = iLengthSignal-iFrameSize+iOverlap;
        iNumberOfFrames = fix( nli/iOverlap  );
    otherwise
        error( 'Undefined method' )
end

% Initialization
vF0   = NaN*ones( 1, iNumberOfFrames );
iT0   = 0;
C     = zeros( 3, iNumberOfFrames );

for iCurrentFrame=1:iNumberOfFrames
    
    % Initial and ending point
    iIni = round( (iCurrentFrame-1)*iOverlap+1 );
    iEnd = iIni-1 + round( iFrameSize );
    
    disp(['Frame ', num2str( iCurrentFrame ), ' of ', num2str( iNumberOfFrames )])
    
    if (iIni>=1) && (iEnd<=iLengthSignal)
        
        % Current frame
        vFrame = vSignal(iIni:iEnd);
        
        switch (sType)
            case 'Boyanov'
                iF0 = PitchBoyanov( vFrame, iFs, iT0 );
            case 'Rabiner_StaticClip'
                % Rabiner -> Static clipping
                iF0 = PitchClip( vFrame, iFs, iT0 );
            case 'Rabiner_LPC'
                % 'Rabiner_LPC'
                iF0 = PitchCorr( vFrame, iFs, iT0 );
            case 'Rabiner_Cepstral'
                % Rabiner -> Cepstral
                iF0 = PitchCeps( vFrame, iFs, iT0);
            case 'Kasuya'
                % Kasuya ->
                iF0 = PitchSegKas( vFrame, iFs );
            otherwise
                error( 'Undefined method' )
        end
        
        if strcmp( sType, 'Kasuya' )
            C(:,iCurrentFrame) = iFs./iF0;
        else
            vF0( iCurrentFrame ) = iF0;
        end
        
        if iF0>0
            % If F0 exists -> Voiced. Then the FrameSize will be updated
            % considering this new value. Otherwise, the segment is
            % unvoiced and the previously calculated FrameSize is used in
            % the next iteration
            iT0    = round(iFs/iF0);
            
            if any( strcmp( sType, {'Boyanov', 'Rabiner_StaticClip'} ) )
                % Updated frameSize
                iFrameSize=iNumberTo*iT0;
            end
        end
    end
end

if strcmp( sType, 'Kasuya' )
    
    % Pitch period, first operating on the median of the first candidates
    % in C
    C1=C(1,:);
    Tmediana=median(C1);
    
    % For each segment, if the first candidate is within 20% of the median,
    % the first candidate is taken as the pitch period of the segment.
    % If not, the same is done with the following candidates. If neither meets the condition,
    % the first is taken as first candidate.
    margen=round(0.2*Tmediana);
    minEntorno=Tmediana-margen;
    maxEntorno=Tmediana+margen;
    
    for n=1:iNumberOfFrames
        if (C(1,n)>minEntorno) && (C(1,n)<maxEntorno)
            vF0(n)=C(1,n);
        elseif (C(2,n)>minEntorno) && (C(2,n)<maxEntorno)
            vF0(n)=C(2,n);
        elseif (C(3,n)>minEntorno) && (C(3,n)<maxEntorno)
            vF0(n)=C(3,n);
        else
            vF0(n)=C(1,n);
        end
    end
    
    vF0( vF0==0 ) = NaN;
end
  
iF0Out = nanmean( vF0 );