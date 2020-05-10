clear variables
clc
close all

addpath( genpath( '../../../src' ) )
addpath( genpath( '../../../External Toolboxes' ) )

sDir = '../../../Audios';

% [vSignalNorm, iFs]  = audioread( fullfile( sDir, '1-a_n.wav' ) );
% vSignalPath  = audioread( fullfile( sDir, '1424-a_n.wav' ) );

[vSignalNorm, iFs]  = audioread( fullfile( sDir, 'asra.wav' ) );
vSignalNorm = normalize( vSignalNorm, 'zscore' );

vSignalPath  = audioread( fullfile( sDir, 'cgra.wav' ) );
vSignalPath = normalize( vSignalPath, 'zscore' );


iFrame  = ceil( 40e-3*iFs );
iSolape = ceil( 0.5*iFrame );
iNumPtos = 100;
sTipo = 'R';

xLab    = 'Frequency (kHz)'; 
yLab    = 'Frequency (kHz)'; 

titulo  = '';
axisIn  = [];
lineW   = 3;
fontS   = 35;

%% Norm
GNE( vSignalNorm, iFs, 1, length( vSignalNorm ), iNumPtos, 50e-3, 500, 100  )

name = 'GNE_Norm'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )

%% Path
GNE( vSignalPath, iFs, 1, length( vSignalPath ), iNumPtos, 50e-3, 500, 100  )

name = 'GNE_Path'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )