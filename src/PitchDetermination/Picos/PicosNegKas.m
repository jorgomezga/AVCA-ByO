function [vAmplitudes, vPosiciones, vPicos]=PicosNegKas(vSignal, iFs, iInicio, iFinal, T)

% Computes the amplitudes P and the positions i_P of the negative peaks, 
% cycle by cycle, of the input voice; it also returns the signal peaks 
% Of the same duration as s
%
% Input parameters:
%   vSignal:    contains the complete speech signal in row or column vector
%   Fs:         the sampling frequency
%   iStart:     first sample of the section to be analyzed
%   iFinal:     last sample of the section to be analyzed
%   T:          the sequence of pitch period values (in samples) calculated
%               on 40ms windows, with 20ms offsets, as they do
%               Kasuya and Feijoo. It is the result of applying the PitchKas function.
%
% Output parameters:
%   P:          amplitudes of the negative peaks, cycle to cycle, of the
%               speech sample
%   i_P:        positions of the negative peak, cycle by cycle, of the voice
%               Peaks of the same duration as vSignal with the peaks found


if nargin<2, iFs=25000; end
if nargin<3, iInicio=1; end
if nargin<4, iFinal=length( vSignal ); end

% Si el vector no es de tipo columna es de suponer que se ha pasado como tipo fila 
if size( vSignal, 2 ) ~= 1, vSignal=vSignal'; end; 

iLongitud=length(T);

% La búsqueda de vPicos se hace a lo largo de todo el tramo de análisis.
% El punto de comienzo es el primer cruce por cero. Desde ese punto, y en un intervalo
% de T(1) muestras, se busca el primer pico negativo. Las siguientes búsquedas se hacen
% desde el anterior pico negativo encontrado, mirando a distancia T(n) muestras, y en un
% entorno 2RT(n). T(n) es el periodo de pitch para el segmento en que estemos situados;
% como hay solape de ventanas consideramos que: segmento(1)=zona_sin_solape_ventana(1),
% segmento(2)=zona_solapada_ventanas (1) y (2), etc. Si un segmento es sordo se salta.

vAmplitudes=[];
vPosiciones=[];
sin_solape=0.02*iFs;

% k es un índice para referenciar a los vectores vAmplitudes (amplitudes de los vPicos negativos), y
% vPosiciones (posiciones de dichos vPicos sobre la señal completa vSignal).
k=1;

% n sirve para indicar el segmento en el que estamos, para utilizar el valor de periodo de 
% pitch T(n) adecuado.
n=1;

R=0.2;

% muestra de referencia para la búsqueda
mueBusq=iInicio;

% comienzo vale 1 si comenzamos la búsqueda desde el primer cruce por cero (esto ocurre si 
% estamos buscando el primer pico, o si acabamos de saltar un segmento sordo), y vale 0 si se
% busca desde el último pico encontrado, a T(n) muestras de distancia en un entorno 2RT(n).
comienzo=1;
terminado=0;

% Si hay segmentos sordos se saltan
while (T(n)==0) && (n<iLongitud),
   n=n+1;
   mueBusq=mueBusq+sin_solape;
end

if (n==iLongitud) && T(n)==0,
   terminado=1;
end

while (mueBusq+T(n)<=iFinal) && (terminado==0),
   
   if comienzo==1,
      % Búsqueda del punto de comienzo:
      if vSignal(mueBusq)==0,
         mueBusq=mueBusq;
      elseif vSignal(mueBusq)>0,
         while vSignal(mueBusq)>0 && (mueBusq<=iFinal),
            mueBusq=mueBusq+1;
         end
      else
         while vSignal(mueBusq)<0 && (mueBusq<=iFinal),
            mueBusq=mueBusq+1;
         end
      end
   
      % Búsqueda del primer pico (o del primero después de un segmento sordo):
      if mueBusq+T(n)<=iFinal,
         [vAmplitudes(k),pos]=min(vSignal(mueBusq:mueBusq+T(n)));
         vPosiciones(k)=pos+(mueBusq-1);
         % El pico encontrado sirve de referencia para la búsqueda del siguiente:
         mueBusq=vPosiciones(k);
         k=k+1;
         comienzo=0;
      end
      
   else
      % El pico encontrado es la muestra de iInicio para la siguiente búsqueda; desde él
      % se mira a T(n) muestras siendo n el número de segmento en que está el pico.
   
      margen=round(R*T(n));
      mueIniBusq=mueBusq+T(n)-margen;
      mueFinBusq=min(mueBusq+T(n)+margen, iFinal);
   
      [vAmplitudes(k),pos]=min(vSignal(mueIniBusq:mueFinBusq));
      vPosiciones(k)=pos+(mueIniBusq-1);
      % El pico encontrado sirve de referencia para la búsqueda del siguiente:
      mueBusq=vPosiciones(k);
      k=k+1;

      n=ceil( (mueBusq-iInicio)/sin_solape );
      if n>iLongitud,
         n=iLongitud;
      end
         
      % Si el segmento es sordo se salta
      while (T(n)==0) && (n<iLongitud),
         comienzo=1;
         n=n+1;
         mueBusq=mueBusq+sin_solape;
      end
         
      if (n==iLongitud) && T(n)==0,
         terminado=1;
      end
      
   end
   
end


% La señal que contiene los vPicos negativos será:
Ns=length( vSignal );
vPicos=zeros( 1,Ns );
if length( vAmplitudes ) > 0,
   vPicos( vPosiciones )=vAmplitudes;
end

if nargout == 0, 
      eje_t=(iInicio:iFinal)/iFs;
      plot(eje_t,vSignal(iInicio:iFinal),'r',eje_t,vPicos(iInicio:iFinal),'g');
      title('Picos negativos');
end; 

return; 

% NOTA: no se ha considerado la posibilidad de que un segmento sea sordo (T(n)=0).


   
   
   
   
      
