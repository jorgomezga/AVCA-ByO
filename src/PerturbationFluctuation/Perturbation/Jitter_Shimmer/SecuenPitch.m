function [vPPS, vPAs, vPitchCont, vAmpCont]=SecuenPitch( vSignal, iFs, iInicio, iFinal, iMetodo)

% Calculates the pitch or amplitude sequences given a vSignal
%
% Input parameters:
% 	vSignal: column vector containing the full speech signal
% 	iFs:	 the sampling frequency
% 	i:	 Start the first sample of the section to be analyzed
% 	iFinal:	 the last sample of the section to be analyzed
% 	iNumPuntos: 	number of HNR values ​​returned
% 	lMetodo:	analysis to be used: 0 (Boyanov), 1 (Kasuya).
%
% Output parameters:
% 	vPPS sequence of pitch periods (cycle to cycle).
% 	vPAS sequence of peak amplitudes of signal cycles. (pitch amplitudes)
% 	vPitchCont sequence "contour of the pitch period" (of the same duration
% 			that the analyzed section; the period in seconds of each cycle
% 			of voice is repeated for all samples of the cycle)
% 	vAmpCont sequence "contour of pitch amplitudes" (of the same duration
% 			that the analyzed section; the amplitude of the peak of each cycle
% 			of voice is repeated for all samples of the cycle)

if nargin < 2, iFs = 50000; end
if nargin < 3, iInicio = 1; end
if nargin < 4, iFinal = length( vSignal ); end
if nargin < 5, iMetodo = 0; end

if iMetodo==0 % (Boyanov)
    % Se obtiene la secuencia de periodos de pitch por el método de Boyanov (ventanas
    % de 3 To, primera ventana de 34 ms, no hay solapamiento entre ventanas).
    T=Pitch( vSignal, iFs, iInicio, iFinal);
    
    % Se eliminan los valores nulos que se asignan al principio y al iFinal si la
    % primera muestra de la señal vSignal o la última muestra se incluyen en el análisis.
    % (El valor 0 se asigna arbitrariamente si la ventana de análisis se sale de
    % la señal de voz. Recuérdese que se toman ventanas centradas en el punto de
    % análisis).
    % T=quita_extremos(T);
    
%     Se obtienen los picos positivos de cada periodo del tramo de señal considerado.
%     (Si se usó PitchKas hay que usar PicosPosKas)
    [vAPos, i_APos, vPicos]=PicosPos( vSignal, iFs, iInicio, iFinal, T);
%     Se obtienen los picos negativos de cada periodo del tramo de señal considerado.
%     (Si se usó PitchKas hay que usar PicosNegKas)
    [vANeg, i_ANeg, vPicos]=PicosNeg( vSignal, iFs, iInicio, iFinal, T);
    
    
%     % Postive peaks
%     [vApos, i_APos, vPicos]=Picos( vSignal, iFs, iInicio, iFinal, vT, 1);
%     
%     % Negative peaks
%     [vAneg, i_ANeg, vPicos]=Picos( vSignal, iFs, iInicio, iFinal, vT, -1);
end

if iMetodo==1 % (kasuya)
    % Se obtiene la secuencia de periodos de pitch por el método de Kasuya (Feijoo), es
    % decir, con ventanas de 40ms y 50% de solape.
    T=PitchKas( vSignal, iFs, iInicio, iFinal);
    
    % Se eliminan los valores nulos que se asignan al principio y al iFinal si la
    % primera muestra de la señal vSignal o la última muestra se incluyen en el análisis.
    % (El valor 0 se asigna arbitrariamente si la ventana de análisis se sale de
    % la señal de voz. Recuérdese que se toman ventanas centradas en el punto de
    % análisis).
    T=QuitaExtremos( T );
    
    % Se obtienen los picos positivos de cada periodo del tramo de señal considerado.
    [vAPos, i_APos, vPicos]=PicosPosKas( vSignal, iFs, iInicio, iFinal, T);
    % Se obtienen los picos negativos de cada periodo del tramo de señal considerado.
    [vANeg, i_ANeg, vPicos]=PicosNegKas( vSignal, iFs, iInicio, iFinal, T);

    
%     [vAPos, i_APos, vPicos]=PicosPosKas( vSignal, iFs, iInicio, iFinal, T, );
%     
%     [P, i_P, picos]=PicosKas( vSignal, iFs, iInicio, iFinal, T, iPosNeg )
end


% Nos quedaremos con los picos que tengan mayor energía
Epos=sum( vAPos )/length( vAPos );
Eneg=abs( sum( vANeg ) )/length( vANeg );

if Epos >= Eneg
    vA=vAPos;
    i_A=i_APos;
else
    vA=abs( vANeg );
    i_A=i_ANeg;
end

iNumPicos=length( vA );     % = length(i_A)

% vA contiene las amplitudes de los picos encontrados, por lo que:
vPAs=vA;

% En i_A están las posiciones, sobre la secuencia de voz vSignal, de los picos encontrados,
% por lo que calculando las distancias entre ellos tendremos la secuencia de periodos
% de pitch en muestras, y dividiendo por iFs en segundos. (Siempre que todo el tramo sea
% sonoro). Además se calculan los contornos de pitch y de amplitudes.
% Si hay segmentos sordos en el tramo de voz considerado, entonces en dichos segmentos
% no se habrán localizado picos, por lo que la distancia entre el último pico antes de
% un segmento sordo y el primero después, no corresponde al valor de un periodo de pitch.
% Consideraremos distancias no válidas si son superiores al máximo periodo admitido, es
% decir 15 ms. En ese caso los contornos toman el valor arbitrario 0.
rPPSmax=0.015;

vPPSTemp=zeros( 1, iNumPicos-1 );

vPitchCont=zeros( 1,iFinal-iInicio+1 );
vAmpCont=zeros( 1, iFinal-iInicio+1 );

k=1;
for i=1:iNumPicos-1,
    vPPSTemp( i )=( i_A(i+1)-i_A(i) )/iFs;
    if vPPSTemp( i ) <= rPPSmax,
        % vPPS contiene sólo los periodos de pitch válidos (menos de 15ms)
        vPPS( k ) = vPPStemp( i );
        k=k+1;
        
        vPitchCont( i_A(i):i_A(i+1)-1 )=vPPSTemp( i );
        vAmpCont( i_A(i):i_A(i+1)-1 )=vPAs( i );
    else
        vPitchCont( i_A(i):i_A(i+1)-1 )=0;
        vAmpCont( i_A(i):i_A(i+1)-1 )=0;
    end
end

% Finalmente se extienden los extremos de los contornos calculados hasta aquí, para
% ocupar todo el intervalo que va de iInicio a iFinal.
vPitchCont( iInicio:i_A(1)-1 )=vPitchCont( i_A(1) );
vPitchCont( i_A(iNumPicos):iFinal )=vPitchCont( i_A(iNumPicos)-1 );

vAmpCont( iInicio:i_A(1)-1 )=vAmpCont( i_A(1) );
vAmpCont( i_A(iNumPicos):iFinal )=vAmpCont( i_A(iNumPicos)-1 );