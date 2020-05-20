function [mParamHist]=histoParamEM(vmMS, bTuned, vMf, vAf)

% [mParamHist]=histoParamEM(vmMS)
% Calculates new parameters starting from Modulation Spectra of an audio
% sequence.
%   INPUTS:
%           - vmMS: Matrix vector containing Modulation Spectra
%   
%   OUTPUTS:
%           - mParamHist: matrix containing in each row histogram-based
%           parameters in the following order:
%                   Cummulative intersection-level (CIL)
%                   Number of points above linear average (PALA)
%                   Ratio of points above linear average (RALA)
%                   RALP25, RALP75, RALP95: Ratio of points above
%                   percentiles (25, 75 and 95)
%


if nargin==1, bTuned=0; end

nudge=1; %histogram bin size. Maximum error when calculating CIL is nudge/2
vEdges1=(-120:nudge:10); % Typical range of Modulation Spectra
vEdges=[-Inf, vEdges1, Inf];

% Indices and variables
iNumHist=length(vEdges);
iNumTramas=size(vmMS,3);
mHist=zeros(iNumTramas,iNumHist);
vPala=zeros(iNumTramas,1);

vRalp25=zeros(iNumTramas,1);
vRalp75=zeros(iNumTramas,1);
vRalp95=zeros(iNumTramas,1);



  %% Percentiles
for i=1:iNumTramas
    
    
    mMS=abs(vmMS(:,:,i));
    
    if bTuned
        mMS25=sintonizaEM(mMS,vMf,vAf,80,0,1800);
        mMS75=sintonizaEM(mMS,vMf,vAf,200,0,2000);
        mMS95=sintonizaEM(mMS,vMf,vAf,200,0,9000);
        mMCIL=sintonizaEM(mMS,vMf,vAf,80,0,2000);
        mMRALA=sintonizaEM(mMS,vMf,vAf,200,800,6000);
       
    else
        mMS25=mMS;
        mMS75=mMS;
        mMS95=mMS;
        mMCIL=mMS;
        mMRALA=mMS;
    end
    
     [iXMS,iYMS,~]=size(mMRALA);
    
    mTemp=20*log10( mMCIL ); %Decibels
    mHist(i,:)=histc(mTemp(:),vEdges); %Histogram
    rMediaLin=mean(mMRALA(:));
    
    

    vRalp25(i)=prctile(mMS25(:),25);
    vRalp75(i)=prctile(mMS75(:),75);
    vRalp95(i)=prctile(mMS95(:),95);
    
    vPala(i)=sum(sum(mMRALA>rMediaLin)); % PALA
        
end

%% CIL 
mCumCreciente=cumsum(mHist,2);
mCumDecreciente=cumsum(fliplr(mHist),2);
mCumDecreciente=fliplr(mCumDecreciente);

vCil=zeros(iNumTramas,1);
for i=1:iNumTramas
    
    [vAux,vYAux]=intersectcurves(vEdges1(1:end-1)+nudge,mCumCreciente(i,2:end-2),...
     mCumDecreciente(i,2:end-2));
    

vCil(i)=vAux(1); %If there are several points, only the first one to intersect is selected.

% figure();
% bar(vEdges1(1:end-1)+nudge,mHist(i,2:end-2),'barWidth',0.5)
% hold on
% plot(vEdges1(1:end-1)+nudge,mCumCreciente(i,2:end-2)', 'Color','m','LineWidth',3)
% title('Modulus')
% plot(vEdges1(1:end-1)+nudge, mCumDecreciente(i,2:end-2)','Color', 'g','LineStyle','--','LineWidth',3)
% plot( [vAux(1) vAux(1)],[0, vYAux(1)], 'k--' )

end

%% RALA & PALA
vRala=vPala./(iXMS*iYMS-vPala);

%Normalizing PALA
vPala=vPala/(iXMS*iYMS);


%% mParamHist

mParamHist=[vCil vPala vRala vRalp25 vRalp75 vRalp95];

