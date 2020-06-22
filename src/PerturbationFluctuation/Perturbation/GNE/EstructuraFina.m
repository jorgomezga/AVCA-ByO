function vEstructuraFina=EstructuraFina( vSenal, iFs, iOrdenLPC, rLongVentana, rDesplazamiento, cWin )

% Calculate the fine structure or the glottic excitation signal (FineStructure)
% of the voice signal (of sampling frequency iFs) by the LPC method
% of iOrdenLPC order autocorrelation, with the cWin window of
% length rLongWindow (in ms) and jump between value windows
% iDisplacement (in ms).
%
% Input parameters:
%   vSenal:     Row vector with the voice signal to which it is going to
%               smooth the spectrum (deconvolve with respect to formants of the vocal tract).
%   iFs:            vSenal sampling rate (Hz)
%   iOrdenLPC:      LPC reverse filter order. Must be less than the length of the window
%   rLongWindow:    Duration of the analysis window
%   rDisplacement:  Duration of the displacement (or jump) between windows
%                   of analysis. Must be less than window length
%   cWin:           String that indicates the type of the analysis window
%
% Output parameters:
%   vFineStructure: ROW vector that stores the result: filtered signal
%                   with the approximation of the inverse of the vocal tract
%                   (fine structure or excitation of the voice signal)

if nargin < 6, cWin='hann'; end
if nargin < 5, rDesplazamiento=10e-3; end
if nargin < 4, rLongVentana=30e-3; end 
if nargin < 3, iOrdenLPC=13; end
if nargin < 2, iFs=10000; end

vEstructuraFina=[];
     
% Normalize the amplitude of the vSignal
vSenal=14000/max(abs(vSenal))*vSenal;
% Remove components around d.c.
vSenal=filtfilt([1 -1],[1 -.99],vSenal);
% do a fixed pre-emphasis
vSenal=filtfilt([1 -0.9], 1, vSenal );
    
iLongVentana=ceil( iFs*rLongVentana );          % In samples
iDesplazamiento=ceil( iFs*rDesplazamiento );    % In samples

iLongFiltroLPC=iOrdenLPC+1;     % Length of the LPC inverse filter

% Overlapping between frames using the overlap-save method
iSolapamiento=iLongVentana-iDesplazamiento; 

if (iOrdenLPC>=iLongVentana)
     error('Error: The order of the LPC filter should be lesser than that of the window ');
end
 
if iSolapamiento < iLongFiltroLPC-1 % Comprueba la validez del parametro iDesplazamiento de entrada
    disp('Error: iDesplazamiento should be <=iLongVentana-(iLongFiltroLPC-1)');
    return;
end

% Change to row vector
if size(vSenal,2)==1, vSenal=vSenal'; end


% OVERLAP-SAVE 
% Add trailing zeros at the beginning for the first vSenal frame to slice
% It is done with the filter, for which the vSenal windows is filled with 
% iDisplacement = N-iLongWindow zeros to have N points output.
% Also with conv, but only FIR -then there would be no need to fill with zeros-
vSenalRell=[zeros(1,iSolapamiento) vSenal]; 

iLongSenalRell=length(vSenalRell);
iNumVentanasRell=ceil((iLongSenalRell-iLongVentana)/iDesplazamiento +1); 
iLongFinalSenalRell=(iNumVentanasRell-1)*iDesplazamiento+iLongVentana;
% Length counting the last incomplete window completed with zeros
% iNumVentanas is the number of windows

% Filled with zeros so that enframe.m considers all the samples
vSenalRell=[vSenalRell zeros(1,iLongFinalSenalRell-iLongSenalRell)]; 

% Windowing. As row vector
vWin=window(cWin,iLongVentana); 
% vWin=hann(iLongVentana);
vWin=vWin'; 

mVentanasLPC=enframe( vSenalRell, vWin, iDesplazamiento ); 
% mLPCWindows: Matrix that stores in each row a windowed frame cWin (Hanning)
% to calculate the LPC coefficients

mVentanas=enframe( vSenalRell, iLongVentana, iDesplazamiento ); 
% mWindows: Array that stores in each row a vSenal frame (rectangular window)

%% LPC calculation and inverse filtering
% We take null initial conditions; You can also start on a later sample of the 
% signal to filter and consider those earlier for initial conditions.
vZi=[]; 
for i=1:iNumVentanasRell
    % Check that the complete window is not filled with zeros, not to divide by zero
    if max(abs(mVentanasLPC(i,:)))==0 
        % Gives the minimum non-null value to the last sample in case of a
        % window filled with zeros
        mVentanasLPC(i,end)=eps;
    end
    % Calculation of the LPC and G coefficients on the Hanning windows; 
    % using the lpc function of MATLAB
    [vCoefFiltroInverso,rGlpc]=lpc( mVentanasLPC(i,:), iOrdenLPC ); 
    vCoefFiltroInverso=real( vCoefFiltroInverso ); 
    
    % The vZf given by filter are not used because there is overlap between windows
    mVentanasFiltradas(i,:)=filter( vCoefFiltroInverso, rGlpc, mVentanas(i,:), vZi ); 
    % Final conditions of this frame, initial conditions for the next
    % iteration
    vZi=filtic( vCoefFiltroInverso, rGlpc, fliplr( mVentanasFiltradas(i,1:iDesplazamiento)), ...
        fliplr(mVentanas(i,1:iDesplazamiento) ) ); 
end

% mFilteredWindows: Matrix with the same structure as windows, with the
% result of the inverse filter of each frame of the resampled voice signal

% Now the fine structure or excitation signal must be reconstructed from each filtered window

% RECONSTRUCTION: Excitation signal
% Overlap-save
% vFineStructure: Signal resulting from reverse filtering: fine structure or
% excitation signal. iLongFinalSenalRell: size
vEstructuraFina=DeFrame( mVentanasFiltradas, iLongVentana, iDesplazamiento );
% Entered zeros are removed to use all vSenal points (truncated)
if ( iLongFinalSenalRell > iLongSenalRell ) 
   vEstructuraFina=vEstructuraFina( 1:iLongSenalRell );
end

% The first overlap samples are removed from the first output frame,
% which are invalid samples (overlap-save)
vEstructuraFina=vEstructuraFina( iSolapamiento+1:end ); 

% Remove components around d.c.
vEstructuraFina=filtfilt([1 -1],[1 -0.99], vEstructuraFina );

% Normalize the amplitude of the analyzing vSignal
vEstructuraFina=14000/max(abs(vEstructuraFina))*vEstructuraFina; 