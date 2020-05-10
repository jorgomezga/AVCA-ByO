function [vParametros, caNombres]=PerturbationFluctuation( vFrame, iFs )

% Calculation of the perturbation and fluctuation parameters
%
% Inputs:
%   vFrame               = Input vFrame
%   iFs                  = Sampling frequency
%
% Outputs:
%   vParametros          = Pitch value of the current frame
%
% Los par�metros contenidos en el vector P son los siguientes:
%     P=[Fo To HNR CHNR NHR VTI SPI NNE jita jitt RAP PPQ sPPQ ShdB shim
%        APQ sAPQ Fftr FTRI Fatr ATRI GNE IB HFEUBI LFEUBI ILFEUBI MN MPN K IK];
% Fo:      frecuencia fundamental media (Hz)
% To:      periodo fundamental medio (s)
% HNR:     Relaci�n arm�nico a ruido seg�n el m�todo de Yumoto (en dB)
% CHNR:    Relaci�n arm�nico a ruido seg�n el m�todo de Guus de Krom (en dB)
% NHR:     Relaci�n ruido a arm�nico en los m�rgenes de frecuencia que usa la Kay, midiendo el ruido en el cepstrum (en tantos por uno)
% VTI:     Indice de turbulencia de la voz de la Kay. La medida se hace utilizando el cepstrum. (En tantos por uno)
% SPI:     Indice de fonaci�n suave de la Kay. La medida se hace utilizando el cepstrum. (En tantos por uno)
% NNE:     Normalized Noise energy
% GNE:     Glottal to Noise Excitation Ratio 
% jita:    jitter absoluto medio (microsegundos)
% jitt:    jitter relativo respecto al periodo fundamental medio (%)
% RAP:     Perturbaci�n media relativa del periodo fundamental (%)
% PPQ:     Cociente de perturbaci�n del periodo fundamental (%)
% sPPQ:    PPQ suavizado. Utiliza un factor de suavizado de 55 ciclos glotales (%)
% ShdB:    shimmer medio (dB)
% shim:    shimmer relativo a la amplitud de pico media (%)
% APQ:     Cociente de perturbaci�n de la amplitud (%)
% sAPQ:    APQ suavizado. Usa un factor de suavizado de 55 ciclos glotales (%)
% Fftr:    frecuencia del tremor de Fo (Hz)
% FTRI:    Indice de intensidad del tremor de Fo (%)
% Fatr:    frecuencia del tremor de amplitud (Hz)
% ATRI:    Indice de intensidad del tremor de amplitud (%)
% IB:      la interferencia del indice de bicoherencia
% HFEUBI:  el valor relativo de energia a alta frec. del indice de bicoherencia unidimensional 
% LFEUBI:  el valor relativo de energia a baja frec. del indice de bicoherencia unidimensional
% ILFEUBI: interferencia del valor de energ�a a baja frec. del indice de bicoherencia unidimensional 
% MN:      la variacion del ruido estimado por medio del modulo del bispectrum
% MPN:     la variacion del ruido estimado pero teniendo en cuenta la variacion de fase
% IK:      la interferencia de kurtosis
% PLI:     pathological likelihood index 

iNumPuntos = 100; 

vFrame     = normalize( vFrame, 'zscore' );
iLongitud  = length( vFrame );

% Pitch period calculation
rFo = PitchBoyanov( vFrame, iFs, 0, 'Temporal' );
rTo = 1/rFo;

%% Additive noise
% Harmonics to noise ratio - Yumoto
[vHNR, rHNR]  = HNRYum( vFrame, iFs, 1, iLongitud, iNumPuntos );

% Glottal-to-noise excitation ratio
[vGNE, rGNE]  = GNE( vFrame, iFs, 1, iLongitud, 1000, 300, iNumPuntos );

% Cespstrum harmonics-to-noise ratio
[vCHNR, rCHNR]     = HNR( vFrame, iFs, 1, iLongitud, iNumPuntos);
[vCHNR2, rCHNR2]   = CHNRs( vFrame, iFs, '' );

% Normalized noise energy
[vNNE, rNNE]  = NNE( vFrame, iFs, 1, iLongitud, iNumPuntos );

%% Modulation noise perturbation features

% For the calculation of the disturbances and tremor some information is needed that
% is obtained from the pitch analysis
[vPPS, vPAS] = ParamPitch( vFrame, iFs, 1, iLongitud, 0 );

% Pitch perturbation measures: jita, jitt, RAP, PPQ, sPPQ.
[rJita, rJitt] = jitter( vPPS );
rRAP = rap( vPPS );
rPPQ = ppq( vPPS );
try
    rSPPQ=sppq( vPPS, 55 );
catch
    warning('rSPPQ not calculated');
    rSPPQ=-1;
end

% Amplitude perturbation measures: shima, shim, APQ, sAPQ.
[rShdB, rShim]=shimmer( vPAS );
try
    rAPQ = apq( vPAS );
catch
    warning('rAPQ not calculated')
    rAPQ=-1;
end

try
rSAPQ = sapq( vPAS, 55 );
catch
    warning('rSAPQ not calculated')
    rSAPQ=-1;
end

%% Fluctuation
% Par�metros de medida del tremor: se usa el procedimiento de la Kay; los �ndices de
% intensidad salen muy distintos a los que da la Kay.
try
[rFTRI, rATRI, rFftr, rFatr]=Tremor( vFrame, iFs );
catch
    warning('Tremor parameters not calculated');
   rFTRI=-1 ;
   rATRI=-1 ;
   rFftr=-1 ;
   rFatr=-1 ;
end


%% TO CHECK
% NHR, VTI y SPI de la Kay (calculadas en el dominio cepstral)
%[vNHR, vVTI, vSPI]=NoiseKAY( vSignal, iFs, 1, iLongitud );
[vNHR, vVTI, vSPI]=KAY( vFrame, iFs, 1, iLongitud);

%[rPLI, rPLIN]=PLI(vSignal, iFs, 1024 );

%[rIB, rHFEUBI, rLFEUBI, rILFEUBI, rMN, rMPN, rIK]=hos( vSignal, iFs, 1, length( vSignal ), iNumPuntos);

rNHR=mean( vNHR );
rVTI=mean( vVTI );
rSPI=mean( vSPI );






% Inicalizamos valores de par�metros 
rIB=0; rHFEUBI=0; rLFEUBI=0; rILFEUBI=0; rMN=0; rMPN=0; rIK=0; rPLIN=0; 
%[rIB, rHFEUBI, rLFEUBI, rILFEUBI, rMN, rMPN, rIK]=hos( vSignal, iFs, 1, length( vSignal ), iNumPuntos);

      
vParametros=[rFo rTo rHNR rCHNR rNHR rVTI rSPI rNNE rGNE rJita rJitt rRAP rPPQ rSPPQ rShdB rShim rAPQ rSAPQ rFftr rFTRI rFatr rATRI rIB rHFEUBI rLFEUBI rILFEUBI rMN rMPN rIK rPLIN]; 


caNombres=['Fo' 'To' 'HNR' 'CHNR' 'NHR' 'VTI' 'SPI' 'NNE' 'GNE' 'jita' 'jitt' 'RAP' 'PPQ' 'sPPQ' 'ShdB' 'shim' 'APQ' 'sAPQ' 'Fftr' 'FTRI' 'Fatr' 'ATRI' 'IB' 'HFEUBI' 'LFEUBI' 'ILFEUBI' 'MN' 'MPN' 'IK' 'PLI' ];


% Normalizamos la amplitud de la se�al a analizar 
% vFrame=14000/max(abs(vFrame))*vFrame;
