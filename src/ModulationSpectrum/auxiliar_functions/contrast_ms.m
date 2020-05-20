function [rContraste]=contrast_ms(mInput,iLongMask)

% function [mContraste]=contrast_ms(mInput)
% Contrast in the mInput image calculated using the Weber-Fechner contrast
% relation, as stated in:
%
% [1] R. Peters and R. Strickland, "Image complexity metrics for automatic
% target recognizers", Autom. Target Recognizer Syst. Technol. Conf., 1990.
%
% It is calculated in a 9 pixels square.
%
% INPUTS:
%       - mInput: image to detect contrast
%       - iLongMask: size of the subgroup of pixels (iLongMask x
%       iLongMask) to use in the contrast filter
%
% OUTPUT:
%       - iContraste: Global contrast value.


if nargin<2, iLongMask=1; end



%Image normalization. If there are negative or 0 values, the image is scaled.
if any(any(mInput<1))
    mInput=mInput+abs(min(min(mInput)))+1;
end

mInput=256*mInput/(max(max(mInput)));



% if iLongMask=1, Pixels are considered as regions. Otherwise, regions are
% processed.
if iLongMask>1
    mH=fspecial('average',iLongMask);
    mAverage=filter2(mH,mInput);
    
    
else
    mAverage=mInput;
end
[iFil,iCol]=size(mAverage);



if iFil==2 %Only two lines in mAverage.
    mcj=zeros(iFil-1,iCol-1);
    for j=2:iCol-1
        
        r1=mAverage(1,j);
        r2=mAverage(2,j);
        r3=mAverage(1,j-1);
        r4=mAverage(1,j+1);
        
        mcj(1,j-1)=(r2-r1)/(r2+r1)+(r3-r1)/(r3+r1)+(r4-r1)/(r4+r1);
    end
    
    
elseif iFil==1
    mcj=zeros(1,iCol-1);
    for j=2:iCol-1
        r1=mAverage(1,j);
        
        r3=mAverage(1,j-1);
        r4=mAverage(1,j+1);
        
        mcj(1,j-1)=(r3-r1)/(r3+r1)+(r4-r1)/(r4+r1);
        
    end
    
else
    mcj=zeros(iFil-1,iCol-1);
    for i=2:iFil-1
        for j=2:iCol-1
            
            r1=mAverage(i,j);
            r2=mAverage(i-1,j);
            r3=mAverage(i,j-1);
            r4=mAverage(i,j+1);
            r5=mAverage(i+1,j);
            
            mcj(i-1,j-1)=(r2-r1)/(r2+r1)+(r3-r1)/(r3+r1)+(r4-r1)/(r4+r1)+(r5-r1)/(r5+r1);
            
        end
        
    end
end

rContraste=sum(sum(mcj))/length(mcj(:));
end






