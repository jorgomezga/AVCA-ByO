function iF0 = StaticClip( vFrame, iFs, iT0 )

% Static clipping 
%
% Inputs:
%   vFrame            = Input vFrame
%   iFs               = Sampling frequency
%   iT0               = Previously calculated pitch period
% Outputs:
%   vClippedFrame     = Clipped frame


% Si el segmento anterior era sonoro, teniendo un periodo de pitch ToAnt~=0, se busca el
% periodo de pitch del segmento actual en un entorno del 30% de ToAnt (siempre que no se
% salga del rango 2ms-15ms). Si no se busca el periodo de pitch entre los valores de 2ms
% y 15ms
if iT0~=0
   iTmin=max( round(0.7*iT0), round(0.002*iFs) );
   iTmax=min( round(1.3*iT0), round(0.015*iFs) );
else
   iTmin=round( 0.002*iFs );
   iTmax=round( 0.015*iFs );
end

UMBRAL_ENER=20; 

iLongitud=length( vFrame );
vMeanFrame=mean( vFrame );

% Whitening to cleanse the input signal and obtain a more clear peak of the
% maximum of the autocorrelation function
% The maximum and minimum signal amplitudes in the first and last thirds of the signal 
% are computed, using a 65% positive clipping level of the minimum value of
% the obtained maxima and a negative clipping level of the maximum value of
% the obtained minima
iMax1=max( vFrame(1:fix(iLongitud/3)) );
iMax2=max( vFrame(fix(2*iLongitud/3):iLongitud) );
rAmax=min( iMax1, iMax2 ); 

iMin1=min( vFrame( 1:fix(iLongitud/3)) );
iMin2=min( vFrame( fix(2*iLongitud/3):iLongitud) );
rAmin=max( iMin1, iMin2 );


if rAmin >= vMeanFrame
   rAmin=min( vFrame );
end

if rAmax <= vMeanFrame
   rAmax=max( vFrame );
end


rCLpos=0.65*( rAmax-vMeanFrame );
rTRmax=rCLpos+vMeanFrame;

rCLneg=0.65*( vMeanFrame-rAmin );
rTRmin=-rCLneg+vMeanFrame;

vSigClip=zeros( 1, iLongitud ); 

for i=1:iLongitud
   vSigClip(i)=vMeanFrame; 
end

for i=1:iLongitud   
   if vFrame(i) > rTRmax
      vSigClip(i)=vFrame(i);
   elseif vFrame(i) < rTRmin
      vSigClip(i)=vFrame(i);      
   end    
end
   
vAutocorr=xcorr( vSigClip );
vAutocorr=vAutocorr( iLongitud+iTmin:2*iLongitud-1 ).^2;

[rMaximo, iMaximo]=max( vAutocorr( 1:iTmax-iTmin ));
   
% Inicializamos el pithc a 0. Como si fuera no sonora 
iF0=0;
if ( rMaximo >UMBRAL_ENER*mean( vAutocorr ) ) 
		
   % Ver un ejemplo sencillo de autocorrelaci�n con periodicidad T
   % Nota: R(tam_vent)=R(0)
   iT=iMaximo+iTmin-1;
   if iT~=0
      iF0=iFs/iT;
   end   
end