function vClippedFrame = DynamicClip( vFrame )

% Dynamic central clipping with threshold adaptation
%
% Inputs:
%   vFrame               = Input vFrame
% Outputs:
%   vClippedFrame        = Clipped frame

iLengthFrame = length( vFrame );
vMeanFrame   = mean( vFrame );

% The maximum and minimum signal amplitudes in the first and last thirds of the signal 
% are computed, using a 65% positive clipping level of the minimum value of
% the obtained maxima and a negative clipping level of the maximum value of
% the obtained minima
max1=max( vFrame( 1:fix( iLengthFrame/3 ) ) );
max2=max( vFrame( fix( 2*iLengthFrame/3 ):end) );
Amax=min( max1, max2 );
if Amax<=vMeanFrame
    Amax=max( vFrame );
end
CLpos=0.65*(Amax-vMeanFrame);
TRmax=CLpos+vMeanFrame;

min1=min( vFrame( 1:fix( iLengthFrame/3 ) ) );
min2=min( vFrame( fix( 2*iLengthFrame/3 ):end ) );
Amin=max( min1, min2 );
if Amin>=vMeanFrame
    Amin=min( vFrame );
end
CLneg=0.65*(vMeanFrame-Amin);
TRmin=-CLneg+vMeanFrame;

% Values to sum or substract to the thresholds in case that
% corrections are needed
ppos=CLpos/5;
pneg=CLneg/5;

% Initialization of variables for the dynamic clipping algorithm.
% ROTmax (or min) is the relationship between positive (or negative) samples
% that are larger than the positive (or negative) threshold with respect to
% the total number of samples (NmaxOK/Ntot).
% ROT1max (or min) and ROT2max (or min) are the inferior and superior
% extremes respectively of the interval that has been proper for
% consideration in the clipping operation
ROT1max=0.2;
ROT2max=0.6;
ROT1min=0.2;
ROT2min=0.6;

clip_valido=0;
clip_max_valido=0;
clip_min_valido=0;

% To avoid oscillations
SentidoPos=0;
SentidoNeg=0;

seg_clip=zeros(1,iLengthFrame);
while clip_valido==0
    for i=1:iLengthFrame
        seg_clip(i)=vMeanFrame;
    end
    
    NmaxOK=0;
    NminOK=0;
    for I=1:iLengthFrame
        if vFrame(I) > TRmax
            seg_clip(I)=vFrame(I);
            NmaxOK=NmaxOK +1;
        elseif vFrame(I) < TRmin
            seg_clip(I)=vFrame(I);
            NminOK=NminOK +1;
        end
    end
    
    if clip_max_valido==0
        ROTmax=NmaxOK/iLengthFrame;
        if ROTmax<ROT1max
            
            if SentidoPos==1
                % To avoid oscillations the clipping is finished
                clip_max_valido=1;
            else
                SentidoPos=-1;
                TRmax=TRmax-ppos;
            end
            
        elseif ROTmax>ROT2max
            
            if SentidoPos==-1                
                % To avoid oscillations the clipping is finished
                clip_max_valido=1;
            else
                SentidoPos=1;
                TRmax=TRmax+ppos;
            end
        else
            clip_max_valido=1;
        end
    end
    
    if clip_min_valido==0        
        ROTmin=NminOK/iLengthFrame;
        
        if ROTmin<ROT1min
            if SentidoNeg==-1
                % To avoid oscillations the clipping is finished
                clip_min_valido=1;
            else
                SentidoNeg=1;
                TRmin=TRmin+pneg;
            end
            
        elseif ROTmin>ROT2min
            if SentidoNeg==1
                % To avoid oscillations the clipping is finished
                clip_min_valido=1;
            else
                SentidoNeg=-1;
                TRmin=TRmin-pneg;
            end
            
        else
            clip_min_valido=1;
        end
    end
    
    if (clip_max_valido==1) && (clip_min_valido==1)
        clip_valido=1;
    end
end

vClippedFrame=seg_clip;