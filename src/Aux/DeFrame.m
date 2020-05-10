function vSenal=DeFrame( mTramas, iLongVentana, iDesplazamiento)

% Reconstructs the vSenal signal from the windows frames created by the
% function deframe.m
%
% Input parameters:
% mFrames: Matrix of iNumWindows x iLongWindow (number of
%           windows x window length) with the window frames 
%	    to be reconstructed
% iLongWindow: Frame length
% iDisplacement: Displacement between successive windows
%
% Output parameters:
% vSenal: Complete signal reconstructed with the window frames 
%
% When deconstructing the signal with enframe.m, the last samples that
% do not enter a new window or frame are discarded. Therefore, this function
% will not rebuild the vSenal signal exactly except in the case iLongitudSenal = iLongFinalSenal,
% To avoid this, you must fill the signal with zeros until making iLongitudSenal = iLongitudFinalSenal;
% when used in combination with calculaVentanas for the calculation of the fine structure of
% the voice signal, structureFina.m, this is taken into account.

vSenal=[];
   
[iNumVentanas, iNumColumnas]=size(mTramas);
% iNumVentanas: Number of windows
% iNumColumnas: Length of the window

if iNumColumnas~=iLongVentana
   error('Error: The length of the window does not coincide with the number of frames of the enframed signal');
end

iLongitudSenal=(iNumVentanas-1)*iDesplazamiento+iLongVentana; 
%iLongitudSenal: Total length of the reconstructed signal

vSenal=zeros(1,iLongitudSenal);

% Reconstruction
vSenal(1,1:iLongVentana)=mTramas(1,:);
for i=2:iNumVentanas
    vSenal(1,((i-2)*iDesplazamiento+iLongVentana+1):((i-1)*iDesplazamiento+iLongVentana))=mTramas(i,iLongVentana-iDesplazamiento+1:end);
end