function vSignalRecortada=QuitaExtremos( vSignal, iUmbral )

% This function removes the 0's (or low amplitude values ) at the beginning or at the end 
% of the signal, vSignal, if any.
% Its utility is to eliminate the null values of the windowed parameters (pitch, energy, ...)
% assigned at the beginning and end of the voice signal (when the analysis
% window is outside the boundaries of voice signal).
%
% Input parameters:
%   vSignal: vector that contains the section of signal to analyze
%
% Output parameters:
%   vSignal Trimmed: vector containing the trimmed signal segment

if nargin<2
    iUmbral=0.02;
end

% Check that the vector is of type column 
if ~isvector(vSignal)
    error( 'El parámetro vSignal no es un vector columna' ); 
elseif size( vSignal, 2 ) ~= 1
    vSignal=vSignal'; 
end 

iLongitud=length( vSignal );

% First sample different of zero
i=1;
while abs(vSignal(i))<=iUmbral
   i=i+1;
   if i==iLongitud
       % This occurs when the signal is all filled with 0's
       break 
   end  
end

% Last sample different of zero
f=iLongitud;
while abs(vSignal(f))<=iUmbral
   f=f-1;
   if f==0       
       % This occurs when the signal is all filled with 0's
       break
   end
end

vSignalRecortada=vSignal(i:f);