function fo=PitchDetecBoy( vACorr,iFs, iVoicedThreshold, iT0 )

% Calculates the pitch using the Boyanov algorithm
%
%   vACorr               = Autocorrelation function of a: 1) Central clipped
%                          frame or 2) Cepstrum liftered frame.
%   iFs                  = Sampling frequency
%   iVoicedThreshold     = Threshold that the maximum of vACorr has to surpass
%                          in order to consider the segment as voiced and
%                          calculate pitch
%   iT0                  = Pitch period of the previous frame
%
% Outputs:
%   vFo                  = Pitch value of the current frame

% Since the input frame is real, vACorr is symetrical, and therefore, only half the
% samples are needed
iLengthSignal  = ( length(vACorr)+1 )/2;
vACorr         = vACorr( iLengthSignal:2*iLengthSignal-1 );

% Normalization
vNormACorr = vACorr/vACorr(1);

% If the previous segment was voiced (pitch period To > 0),
% the current segment pitch period is looked for, around the 30% of T0
% (as long as the range does not exceed 2 ms - 15 ms). If not, the
% pitch period is looked in the range between the values of 2 ms and 15 ms.
if iT0>0
   Tmin = max( round(0.7*iT0), round(0.002*iFs) );
   Tmax = min( round(1.3*iT0), round(0.015*iFs) );
else
   Tmin = round(0.002*iFs);
   Tmax = round(0.015*iFs);
end

% Sanity check
Tmin = max( 1, Tmin );
Tmax = min( length(vACorr), Tmax );

% First maximum of the normalized autocorrelation
[Rmax,i_max]=max(vNormACorr(Tmin:Tmax));

% If larger than the iVoicedThreshold -> the frame is voiced.
if Rmax>iVoicedThreshold
    % The position of the first maximum is i_max+(Tmin-1). 
    % The values ds are the distances between the maxima of the
    % autocorrelation function
   maximos(1)=i_max+(Tmin-1);
   ds(1)=maximos(1)-1;
   
   % q is a threshold to fix the ATR Threshold for the search of the
   % maximum values
   q=0.6;
   ATR=q*Rmax;
   
   % Maxima values of vNormACorr which are larger than ATR and are at least
   % to Tmin samples of the previous samples are searched
   i=1;
   fin_busq=0;
   while (fin_busq==0)
      mueIni=Tmin+maximos(i);
      if mueIni>iLengthSignal
         fin_busq=1;
      else
         [Rm,i_m]=max(vNormACorr(mueIni:iLengthSignal));
         if Rm<ATR
            fin_busq=1;
         else
            i=i+1;
            maximos(i)=i_m+(mueIni-1);
            ds(i)=maximos(i)-maximos(i-1);
         end
      end
   end
      
   N_max=length(maximos);
   if N_max==1
       % If only the first maximum is larger than the ATR value, then this is
       % taken as the pitch period in samples.
      To=ds(1);
      if To~=0
         fo=iFs/To;
      end
      
   else
      % If several maxima values are detected, then the distance between
      % them is computed and then it is checked if their values are
      % approximately equal. Otherwise f0 is not computed and f0=0
      i=2;
      distOK=1;
      while (i<=N_max) && (distOK==1)
         if ds(i)-ds(i-1)>0.4*ds(i)
            distOK=0;
         end
         i=i+1;
      end
      
      if distOK==1
         To=ds(1);
         if To~=0
            fo=iFs/To;
         end
      end
   end
   
else    
    % Otherwise the segment is unvoiced and f0 = 0 
    fo=0;
end