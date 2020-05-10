clear variables
close all

% addpath( genpath( '/home/jorge/Documentos/Repositorio' ) )
% 
% sDirNorm = '/mnt/Data/BasesDatos/Saarbruecken/DatabaseSeg_Particion_Janaina_16_a_69_rev3/A/Condicion/Norm';
% sDirPath = '/mnt/Data/BasesDatos/Saarbruecken/DatabaseSeg_Particion_Janaina_16_a_69_rev3/A/Condicion/Path';

addpath( genpath( '../../../src' ) )
addpath( genpath( '../../../External Toolboxes' ) )
addpath( '../' )

sDir = '../../../Audios';

% [vSignalNorm, iFs]  = audioread( fullfile( sDir, '1-a_n.wav' ) );
% vSignalPath  = audioread( fullfile( sDir, '1424-a_n.wav' ) );

[vSignalNorm, iFs]  = audioread( fullfile( sDir, 'asra.wav' ) );
vSignalNorm = normalize( vSignalNorm, 'zscore' );

vSignalPath  = audioread( fullfile( sDir, 'cgra.wav' ) );
vSignalPath = normalize( vSignalPath, 'zscore' );


iFrame    = ceil( 40e-3*iFs ); 
rSolape   = 0.5;
iSolape   = floor( (1 - rSolape)*iFrame );

bSmooth   = 1;
sNormaliza     = 'line';

idx     = 10;
xLab    = 'Quefrency (ms)'; 
yLab    = 'Gamnitude (dB)'; 
titulo  = '';
axisIn  = [100 400 7 13];
lineW   = 3;
fontS   = 30;

%% Path
mSignalPath  = enframe( vSignalPath, hamming( iFrame ), iSolape );
mParametrosCPP = ( cepstral_peak_prominence( mSignalPath', iFs, bSmooth, sNormaliza, idx ) )';

name = 'CPP_Path'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )

% Modify cepstral_peak_prominence (cpp.m in the COVAREP toolbox) to return CepsLim and quefSeq
% and then change idx to the frame you're interested to paint
%
% The cepstral_peak_prominence used here has been modified to accept an
% extra parameter idx, and do the whole computation of mSignal (framing) in
% parallel,
% The graphic presented in the paper should be reproducible using the following code
%
% idx = 10;
% p = polyfit( quefSeq, CepsLim(:,idx), 1 );
% Y = polyval( p, CepsLim(:,idx));
% 
% figure
% plot(quefSeq,CepsLim(:,idx),'k', quefSeq,Y, 'r' );

%% Norm
mSignalNorm  = enframe( vSignalNorm, hamming( iFrame ), iSolape );

mParametrosCPP = ( cepstral_peak_prominence( mSignalNorm', iFs, bSmooth, sNormaliza, idx ) )';
name = 'CPP_Norm'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )