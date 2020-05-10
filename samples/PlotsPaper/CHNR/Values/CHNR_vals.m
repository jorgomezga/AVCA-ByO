clear variables
clc
close all

addpath( genpath( '../../../../src' ) )
addpath( genpath( '../../../../External Toolboxes' ) )

sDir = '../../../../Audios';

% [vSignalNorm, iFs]  = audioread( fullfile( sDir, '1-a_n.wav' ) );
% vSignalPath  = audioread( fullfile( sDir, '1424-a_n.wav' ) );

[vSignalNorm, iFs]  = audioread( fullfile( sDir, 'asra.wav' ) );
vSignalNorm = normalize( vSignalNorm, 'zscore' );

vSignalPath  = audioread( fullfile( sDir, 'cgra.wav' ) );
vSignalPath = normalize( vSignalPath, 'zscore' );

iFrame  = ceil( 40e-3*iFs );
iSolape = ceil( 0.5*iFrame );
sTipo = 'R';

xLab    = 'Time (s)'; 
yLab    = 'CHNR (dB)'; 
titulo  = '';
axisIn  = [0 1 0 40];
lineW   = 3;
fontS   = 35;

%% Norm
mSignal = enframe( vSignalNorm, hamming( iFrame ), iSolape );

for j=1:size( mSignal, 1 )
    vFrame = mSignal(j,:);
    HNR_Norm(j) = HNRi( vFrame, iFs );
end

xVal = linspace( 0, length(vSignalNorm)/iFs-iFrame/iFs, length(HNR_Norm) );
stem( xVal, HNR_Norm )

name = 'HNR_Norm'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )

%% Path
mSignal = enframe( vSignalPath, hamming( iFrame ), iSolape );

for j=1:size( mSignal, 1 )
    vFrame = mSignal(j,:);
    HNR_Path(j) = HNRi( vFrame, iFs );
end

xVal = linspace( 0, length(vSignalPath)/iFs-iFrame/iFs, length(HNR_Path) );
stem( xVal, HNR_Path )

name = 'HNR_Path'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )