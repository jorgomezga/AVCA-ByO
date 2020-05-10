clear variables
close all
clc

addpath( genpath( '../../../../src' ) )
addpath( genpath( '../../../../External Toolboxes' ) )

sDir = '../../../../Audios';

% vSignalNorm  = audioread( fullfile( sDir, '1-a_n.wav' ) );
% vSignalPath  = audioread( fullfile( sDir, '1424-a_n.wav' ) );
[vSignalNorm, iFs]  = audioread( fullfile( sDir, 'asra.wav' ) );
vSignalNorm = normalize( vSignalNorm, 'zscore' );

vSignalPath  = audioread( fullfile( sDir, 'cgra.wav' ) );
vSignalPath = normalize( vSignalPath, 'zscore' );


idx       = 10;
iFrame    = ceil( 40e-3*iFs ); 
rSolape   = 0.5;
iSolape   = floor( (1 - rSolape)*iFrame );


xLab    = '$Log_{10}(L)$';
yLab    = '$\log_{10}{g_{k}^{abs}}$';
titulo  = '';
lineW   = 3;
fontS   = 35;
axisIn  = [0.7 1.4 -0.7 -0.2];

%% Norm
mSignal  = enframe( vSignalNorm, hamming( iFrame ), iSolape );
vFrame = mSignal(idx,:);

iHurst      = hurst_estimate( vFrame, 'absval', 1, 1 );
str = ['$h^e=$', num2str( iHurst, 2 )];

name = 'Hurst_Norm'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS, 1, str )

%% Path
mSignal  = enframe( vSignalPath, hamming( iFrame ), iSolape );
vFrame = mSignal(idx,:);

iHurst      = hurst_estimate( vFrame, 'absval', 1, 1 );
str = ['$h^e=$', num2str( iHurst, 2 )];
% TextLocation( str,'Interpreter','latex', 'Location', 'Best' )

name = 'Hurst_Path'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS, 1, str )