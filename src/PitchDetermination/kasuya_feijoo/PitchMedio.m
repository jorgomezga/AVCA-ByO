function [rFo, vF0]=PitchMedio( vSignal, iFs, iInicio, iFinal )

% Calculates the average pitch in Hz of the input voice 
% The average pitch period is calculated for the voice section under analysis
% using the Kasuya procedure (40ms windows with 50% overlap)
%
% Input parameters
%   vSignal:    vector column containing the complete speech signal
%   iFs:        the sampling frequency
%   i:          Starting sample of the section to be analyzed
%   iFinal:     The last sample of the section to be analyzed
%
% Output parameters
%   rFo:        average pitch calculated over 40ms segments,
%               with displacements of 20ms, as Kasuya and Feijoo do.

if nargin<2, error('Not enough input parameters!'); end;
if nargin<3, iInicio=1; end
if nargin<4, iFinal=length( vSignal) ; end

rFo=0; 

% Si el vector no es de tipo columna es de suponer que se ha pasado como tipo fila 
if size( vSignal, 2 ) ~= 1, vSignal=vSignal'; end; 

vTos=PitchKas( vSignal, iFs, iInicio, iFinal );

% Eliminamos los valores nulos correspondientes a los segmentos sordos
rPromedio=0;
k=0;
for i=1:length(vTos)
   if vTos(i)~=0,
      rPromedio=rPromedio+vTos(i);
      k=k+1;
   end
end

rTmed=0; 
if k~=0, rTmed=rPromedio/k; end; 
if rTmed~=0, rFo=iFs/rTmed; end; 

vF0 = iFs./vTos;