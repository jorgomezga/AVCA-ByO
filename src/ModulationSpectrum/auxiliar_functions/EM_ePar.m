function [ mFeatures, caFeatureNames ] = EM_ePar( ePar, vSignal, sTipo, iSubMod, eOptions, iVerbosity)

% Generates the MS features as detailed in refs [1] and [2] using the MS
% already calculated.
% INPUTS:
%       ePar - Structure containing modulation spectrum in the frames of a
%       signal (vSignal) calculated with the function EM.m
%       vSignal - vector containing the input signal (mono)
%       iFs - The sampling rate of vSignal, in Hz.
%       sTipo    any sensible combination of the following:
%               'R'  rectangular window in time domain
%               'N'  Hanning window in time domain
%               'M'  Hamming window in time domain (default)
%               't'  triangular shaped filters in mel domain (default)
%               'n'  hanning shaped filters in mel domain
%               'm'  hamming shaped filters in mel domain
%               'C' calculate Modulation Spectra Centroids
%               'Y' calculate Modulation Spectra Dynamic Range for every
%                   Modulation Frequency Band
%               'f' normalize Centroids respect to pitch
%               'F' normalize Centroids respect to the 0Hz centroid
%               'O' calculate "Low Modulation Ratio" in Modulation Spectra.
%               'P' calculate Contrast, and homogeneity in Modulation Spectra modulus
%                   and phase
%               'A' Calculate Modulation Spectra Histogram features.
%        iSubMod - Number of centroids
%        eOptions - Structure containing other auxiliar and secondary
%        inputs, such as:
%               - iFmodMax - maximum modulation frequency
%               - vFlmr - max frequency used in LMR calculation
%               - EMiSubbands - desired number of acoustic subbands
%               - EMsSpecopt - desired number of modulation bands
%               - EMvFlmr - High bound of the band used to calculate LMR
%               (25 Hz by default)
%   	 iVerbosity - if iVerbosity==1, the function displays texts related
%   	 to the calculation process in screen
%
% OUTPUTS:
%       mFeatures - Matrix containing the MS features. Each row will
%       be referred to a frame and, depending on the values of sTipo,
%       will contain the features:
%       [Centroids, Dynamic range, LMR, Modulus contrast, Modulus homogeneity, 
%           Phase contrast, Phase homogeneity, CIL, PALA, RALA, RALP25, 
%           RALP75, Ralp95]
%       caFeatureNames - Includes the names of the features calculated. Has
%       the same length as the number of columns of mFeatures. 
%
% REFERENCES:
% [1] Moro-Velázquez, L., Gómez-García, J. A., Godino-Llorente, J. I., & 
% Andrade-Miranda, G. (2015). Modulation spectra morphological parameters: 
% a new method to assess voice pathologies according to the grbas scale. 
% BioMed research international, 2015.
%
% [2] Moro-Velázquez, L., Gómez-García, J. A., & Godino-Llorente, J. I. 
% (2016). Voice pathology detection using modulation spectrum-optimized 
% metrics. Frontiers in bioengineering and biotechnology, 4, 1.

if nargin<1, error('MS structure missing' ); end
if nargin<2, sTipo='M'; end
if nargin<3, iSubMod=26; end
if nargin<4, eOptions=[]; end
if nargin<5, iVerbosity=0; end

caFeatureNames=[];

iFs=ePar.data.fs;
iFrame=ePar.data.origlen;

if ~isempty(eOptions) && isfield( eOptions, 'iShift' )
    iShift=eOptions.iShift;
else
    iShift=ceil(iFrame/2);
end


if isfield (  eOptions, 'EMvFlmr' ), vFlmr=eOptions.EMvFlmr;
else, vFlmr= 25; end

iNumtramas=size(ePar.vmMS,3); % Number of frames


% Audio Source Framing plus windowing
[mSignal]=enframing_avca( vSignal, iFrame, iShift, sTipo );

if size(mSignal,1)~=iNumtramas, error('Audio and parameterization frames differ'); end

%Modulation Spectra Source

vmMS=ePar.vmMS;
vMf=ePar.vMf;
vAf=ePar.vAf;


% If eOptions.EMbTuned == 1, MS features will be calculated considering
% configuration defined as optimal in [2]
if isfield (  eOptions, 'EMbTuned' ), bTuned=eOptions.EMbTuned;
else, bTuned=0; end

if ~bTuned
    if isfield (  eOptions, 'EMiFmodMax' ), iFmodMax=eOptions.EMiFmodMax;
    else, iFmodMax= vMf(end); end
    
    if isfield (  eOptions, 'EMiFaMin' ), iFaMin=eOptions.EMiFaMin;
    else, iFaMin= vAf(1); end
    
    if isfield (  eOptions, 'EMiFaMax' ), iFaMax=eOptions.EMiFaMax;
    else, iFaMax= vAf(end); end
    
    
    if iFaMax>iFs/2, iFaMax=vAf(end);
        % Modulation frequency limit
        if iFmodMax<vMf(end)
            sSpecopt=max(find(vMf<=iFmodMax));
            vmMS=vmMS(:,1:sSpecopt,:);
        end
        
        %Acoustic frequency filtering
        [~,iAcmin]=min((abs(vAf-iFaMin)));
        [~,iAcmax]=min((abs(vAf-iFaMax)));
        
        if iAcmin==1, iAcmin=0;end
        
        vmMS=vmMS(iAcmin+1:iAcmax,:,:);
        vAf=vAf(iAcmin+1:iAcmax);
        
    end
    
end

if isfield ( eOptions, 'EMsSpecopt' ), sSpecopt=eOptions.EMsSpecopt;
else, sSpecopt= size(vmMS,2); end

clear('ePar');

% mParametros will contain the output features
mFeatures=[];

% Centroids Calculation
if any(sTipo=='C')|| any(sTipo=='Y')
    %Modulation frequency band reduction
    red_vmMS = band_reduction(vmMS,iSubMod*2);
    
    % Only positive part of spectrum is used
    red_vmMS = abs( red_vmMS(:,1:iSubMod,: ) );
    
    %Centroid calculation
    if any(sTipo=='C') && any(sTipo=='Y')
        mFeatures=zeros(iNumtramas,iSubMod*2);
    else
        mFeatures=zeros(iNumtramas,iSubMod);
    end
    
    if any(sTipo=='C')
        for i=1:iNumtramas
            red_mMSCaux=red_vmMS(:,:,i);
            mFeatures(i,1:iSubMod)=(vAf*red_mMSCaux)./sum(red_mMSCaux);
        end
        
        caFeatureNames=cell(1,iSubMod);
        caFeatureNames(:)={'Centroid'};
        %-----------------
        % Normalization with respect to fundamental frequency
        
        if any(sTipo=='f')
            if iVerbosity==1
                disp( 'Centroids normalization respect to fundamental frequency.');
            end
            vTo=PitchKas(vSignal, iFs); % pitch Period calculated every 40 ms
            vTo=vTo(find(vTo)); % Periods with 0 samples are deleted
            iFo=mean(iFs./vTo);
            mFeatures(:,1:iSubMod)=mFeatures(:,1:iSubMod)./iFo;
        end
        
        %-----------------
        % Normalization with respect to the 0Hz centroid.
        if any(sTipo=='F')
            if iVerbosity==1
                disp( 'Centroids normalization respect to 0 Hz centroid');
                for i=1:iNumtramas
                    rNormali=mFeatures(i,1);
                    %First centroid is not normalized due to  value
                    %always will be 1.
                    mFeatures(i,2:iSubMod)=mFeatures(i,2:iSubMod)./rNormali;
                end
            end
        end
    end
    
    % Dynamic range calculation
    if any(sTipo=='Y')  
        red_vmMS=20*log10(red_vmMS);
        for i=1:iNumtramas
            mFeatures(i,end-iSubMod+1:end)=dynamic_range_MS(red_vmMS(:,:,i));
        end
        clear red_vmMS
        auxNames=cell(1,iSubMod);
        auxNames(:)={'Dynamic Range'};
        caFeatureNames=[caFeatureNames, auxNames];
    end
    
end

% Low Modulation Ratio calculation

if any(sTipo=='O')
    
    mLMR=zeros(iNumtramas,length(vFlmr));
    
    
    for i=1:iNumtramas
        %Modulation Spectra energy
        mTemp = 20*log10( fftshift( abs( vmMS(:,:,i) ), 2 ) );
        % pitch calculation for each frame
        vTo=PitchKas(mSignal(i,:), iFs); % pitch Period calculated every 40 ms
        vTo=vTo(find(vTo)); % Periods with 0 samples are deleted
        iFo=iFs/mean(vTo);
        
        for j=1:length(vFlmr)
            %iBandaA contains the acoustic band index in which fundamental frequency is within the MS matrix
            iBandaA=find(iFo<vAf,1)-1;
            
            %iBandaM contains the index of the last modulation frequency band below 25 Hz within the MS matrix
            
            iBandaM=find(vFlmr(j)<vMf,1)-1;
            if iBandaM>sSpecopt/2
                iBandaM=floor(sSpecopt/2);
            end
            
            %If fundamental frequency cannot be calculated, the first acoustic band is
            %employed
            if isempty(iBandaA) || iBandaA==0
                warning(['Fundamental frequency not calculated for frame ' num2str(i) '. Using first acoustic band instead'] )
                iBandaA=1;
            end
            
            if mod(sSpecopt,2)   %Odd number of fft lines:
                %Energy central band
                iEnergia0=10*log10(10.^(mTemp(iBandaA,floor(sSpecopt/2))/10)+10.^(mTemp(iBandaA,ceil(sSpecopt/2))/10)+10.^(mTemp(iBandaA,ceil(sSpecopt/2)+1)/10));
                %Energy up to 25Hz band
                iEnergia25=10*log10(sum(10.^(mTemp(iBandaA,(ceil(sSpecopt/2)-iBandaM+1):(ceil(sSpecopt/2)+iBandaM)))));
                
                mLMR(i,j)=iEnergia25-iEnergia0;
                
            else   %Even number of fft lines:
                %Energy central band
                iEnergia0=10*log10(10.^(mTemp(iBandaA,sSpecopt/2)/10)+10.^(mTemp(iBandaA,sSpecopt/2+1)/10));
                %Energy up to 25Hz band
                iEnergia25=10*log10(sum(10.^(mTemp(iBandaA,(sSpecopt/2-iBandaM+1):(sSpecopt/2+iBandaM))./10)));
                
                mLMR(i,j)=iEnergia25-iEnergia0;
                
            end
        end
    end
    
    mFeatures=[mFeatures mLMR];
    
    caFeatureNames=[caFeatureNames, 'LMR'];

end

clear('mSignal');

%Contrast, and homogeneity in Modulation Spectra modulus and phase
if any(sTipo=='P')
    
    vmPhase=angle(vmMS);
    vmPhase(:,:,1)=unwrap(vmPhase(:,:,1),[],2); %The first frame will be the seed.
    vmPhase=unwrap(vmPhase,[],3);
    vmPhase=circshift(vmPhase,[0,floor(sSpecopt/2),0]);
    % vmPhase=DerivatePh( vmPhase ); % first derivative can be used instead of phase itself.
     
    %  iLongMaskHom defines the the filter size of the homogeneity mask
    %  (iLongMaskHom x iLongMaskHom) employed in the convolution to 
    %  calculate MS homogeneity
    
    if ~bTuned
        if isfield (eOptions, 'EMiLongMaskHom' ), iLongMaskHom=eOptions.EMiLongMaskHom;
        else, iLongMaskHom=6; end
        
    else
        iLongMaskHom=6;
    end
    
    %  iLongMaskCon defines the the filter size of the contrast mask
    %  (iLongMaskCon x iLongMaskCon) employed in the convolution to 
    %  calculate MS contrast
    if isfield (eOptions, 'EMiLongMaskCon' ), iLongMaskCon=eOptions.EMiLongMaskCon;
    else, iLongMaskCon=6; end
    
    
    
    for i=1:iNumtramas
        %Modulus
        mTempM = 20*log10( fftshift( abs( vmMS(:,:,i) ), 2 ) );
        %Phase
        mTempP = vmPhase(:,:,i);
        mParamImagen(i,1)=contrast_ms(mTempM,iLongMaskCon);
        if bTuned
            mMSH=sintonizaEM(vmMS(:,:,i),vMf,vAf,80,200,9000);
            
            mTempM = 20*log10( fftshift( abs( mMSH ), 2 ) );
            mParamImagen(i,2)=homogen_ms(mTempM,iLongMaskHom);
            
        else
            mParamImagen(i,2)=homogen_ms(mTempM,iLongMaskHom);
        end
        mParamImagen(i,3)=contrast_ms(mTempP,iLongMaskCon);
        mParamImagen(i,4)=homogen_ms(mTempP,iLongMaskHom);
    end
    
    mFeatures=[mFeatures mParamImagen];
    caFeatureNames=[caFeatureNames, 'MSW Module', 'MSH Module', 'MSW Phase', 'MSH Phase' ];
    clear('vmPhase');
end

% Histogram 6 features (CIL, PALA, RALA, RALP25, RALP75, RALP95)
if any(sTipo=='A')
    
    [mParamHist]=histoParamEM(vmMS, bTuned, vMf, vAf);
    mFeatures=[mFeatures mParamHist];
    caFeatureNames=[caFeatureNames, 'CIL', 'PALA', 'RALA', 'RALP25', 'RALP75', 'RALP95' ];
    
end

clear('vmMS');

end

