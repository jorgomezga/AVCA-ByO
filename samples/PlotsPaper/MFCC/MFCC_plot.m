clear variables
clc
close all

addpath( genpath( '../../../src' ) )
addpath( genpath( '../../../External Toolboxes' ) )
addpath( '../' )

sDir = '../../../Audios';

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.03], [0.05 0.05], [0.05 0.05]);

% [vSignalNorm, iFs]  = audioread( fullfile( sDir, '1-a_n.wav' ) );
% vSignalPath  = audioread( fullfile( sDir, '1424-a_n.wav' ) );

[vSignalNorm, iFs]  = audioread( fullfile( sDir, 'asra.wav' ) );
vSignalNorm = normalize( vSignalNorm, 'zscore' );

vSignalPath  = audioread( fullfile( sDir, 'cgra.wav' ) );
vSignalPath = normalize( vSignalPath, 'zscore' );
    
%% Parameters
iFrame   = ceil( 40e-3*iFs );
iSolape  = ceil( 0.5*iFrame );
sTipo    = [];
iNumCoef = 14;

xLab    = ''; 
yLab    = ''; 

titulo  = '';
axisIn  = [];
lineW   = 3;
fontS   = 25;

%% Norm
vSignalNorm = vSignalNorm(1:25000);
MFCCs( vSignalNorm, iFs, sTipo, iNumCoef, iFrame, iSolape, [], 0 )
    
name = 'MFCC_Norm'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )

%% Path
vSignalPath = vSignalPath(1:25000);
MFCCs( vSignalPath, iFs, sTipo, iNumCoef, iFrame, iSolape, [], 0 )
    
name = 'MFCC_Path'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )