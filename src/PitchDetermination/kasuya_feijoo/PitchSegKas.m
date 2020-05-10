function vT0 = PitchSegKas( vFrame, iFs )

% Calculates the pitch of a given frame using the method described by
% Kasuya in:
%       Kasuya, H., Ogawa, S., Kikuchi, Y., & Ebihara, S. (1986). An 
%       acoustic analysis of pathological voice and its application to the 
%       evaluation of laryngeal pathology. Speech communication, 5(2), 171-181.
%
% Inputs:
%   vFrame               = Input vFrame
%   iFs                  = Sampling frequency
% Outputs:
%   vFo                  = Column vector of size 1,2 or 3, considering the
%                          possible pitch candidates of the frame according
%                          to the algorithm of Kasuya (candidates are ordered)

% Ensure that the vFrame is a column vector
if isrow( vFrame )
    vFrame = vFrame';
end

% Voiced threshold based on xCorr
VOICED_THRESH=0.3;
% Silence/Non-silence decision making based on thresholding
bSound = SilentDetectorThreshold( vFrame );

% If the segment is silent, by convention, three candidates are initialized to zero
% and the operation is finished
if bSound == false
    vT0(1:3)=0;
else
    % Pitch period should be around 2ms - 15ms
    Tmin=round(0.002*iFs);
    Tmax=round(0.015*iFs);
    
    iLengthFrame = length(vFrame);
    
    iMeanFrame = mean(vFrame);
    
    % Whitening to cleanse the input signal and obtain a more clear peak of the
    % maximum of the autocorrelation function
    % The maximum and minimum signal amplitudes in the first and last thirds of the signal
    % are computed, using a 65% positive clipping level of the minimum value of
    % the obtained maxima and a negative clipping level of the maximum value of
    % the obtained minima
    max1=max(vFrame(1:fix(iLengthFrame/3)));
    max2=max(vFrame(fix(2*iLengthFrame/3):iLengthFrame));
    Amax=min(max1,max2);
    if Amax<=iMeanFrame
        Amax=max(vFrame);
    end
    CLpos=0.65*(Amax-iMeanFrame);
    TRmax=CLpos+iMeanFrame;
    
    min1=min(vFrame(1:fix(iLengthFrame/3)));
    min2=min(vFrame(fix(2*iLengthFrame/3):iLengthFrame));
    Amin=max(min1,min2);
    if Amin>=iMeanFrame
        Amin=min(vFrame);
    end
    CLneg=0.65*(iMeanFrame-Amin);
    TRmin=-CLneg+iMeanFrame;
    
    seg_clip = zeros(1, iLengthFrame);
    for i=1:iLengthFrame
        seg_clip(i)=iMeanFrame;
    end
    
    for I=1:iLengthFrame
        if vFrame(I)>TRmax
            seg_clip(I)=vFrame(I);
        elseif vFrame(I)<TRmin
            seg_clip(I)=vFrame(I);
        end
    end
    
    % Autocorrelation and elimination of repeated information
    R=xcorr(seg_clip);
    R=R(iLengthFrame:2*iLengthFrame-1);
    Rn=R/R(1);
    % Max of ACF
    [m,i_m]=max(Rn(Tmin:Tmax));
    
    if m>VOICED_THRESH
        % If larger than VOICE_THRESH -> First candidate
        vT0(1)=i_m+(Tmin-1)-1;
        
        % Checks for peaks within 20% of 2k1 and 0.5k1, which exceed
        % threshold. These peaks will be the second and third candidates
        % ordered according to their amplitude.
        % If the conditions are not met, the value T = 0 is taken.
        % At the end, vector T contains the candidates to be the fundamental
        % period of the ordered segment.
        Cand2OK=0;
        Cand3OK=0;
        
        % 0.5k1:
        margen=0.2*(0.5*vT0(1));
        mueMin=round(0.5*vT0(1)-margen);
        mueMax=round(0.5*vT0(1)+margen);
        [m,i_m]=max(Rn(mueMin:mueMax));
        
        if m>VOICED_THRESH
            vT0(2)=i_m+(mueMin-1)-1;
            Cand2OK=1;
        else
            vT0(2)=0;
        end
        
        % 2k1:
        if 2*vT0(1)<=iLengthFrame
            margen=0.2*(2*vT0(1));
            mueMin=round(2*vT0(1)-margen);
            mueMax=min(round(2*vT0(1)+margen),iLengthFrame);
            [m2,i_m]=max(Rn(mueMin:mueMax));
            
            if m2>VOICED_THRESH
                vT0(3)=i_m+(mueMin-1)-1;
                Cand3OK=1;
            else
                vT0(3)=0;
            end
        else
            vT0(3)=0;
        end
        % If there is no valid peak at 0.5k1 but there is one at 2k1, then the peak at 2k1
        % becomes the second candidate.
        % If there were two valid peaks in 0.5k1 and 2k1, they are ordered in decreasing order of amplitude.
        % (The one with the widest range is the second candidate, and the other is the third candidate).
        if (Cand3OK==1 && Cand2OK==0) || m2>m
            Temp=vT0(3);
            vT0(3)=vT0(2);
            vT0(2)=Temp;
        end
       
    else
        % If the candidate is unvoiced:
        vT0(1:3)=0;
    end
    
end

% T is returned as column vector
vT0=vT0';