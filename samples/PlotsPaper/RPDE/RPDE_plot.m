clear variables
close all
clc

addpath( genpath( '../../../src' ) )
addpath( genpath( '../../../External Toolboxes' ) )

sDir = '../../../Audios';

% vSignalNorm  = audioread( fullfile( sDir, '1-a_n.wav' ) );
% vSignalPath  = audioread( fullfile( sDir, '1424-a_n.wav' ) );

[vSignalNorm, iFs]  = audioread( fullfile( sDir, 'asra.wav' ) );
vSignalNorm = normalize( vSignalNorm, 'range', [0 0.9] );

vSignalPath  = audioread( fullfile( sDir, 'cgra.wav' ) );
vSignalPath = normalize( vSignalPath, 'range', [0 0.9] );


idx       = 3;
iFrame    = ceil( 40e-3*iFs ); 
rSolape   = 0.5;
iSolape   = floor( (1 - rSolape)*iFrame );

m   = 4;
tau = 10;

xLab    = '$p(\Delta t_r)$';
yLab    = '$\Delta t_r$';
titulo  = '';
lineW   = 3;
fontS   = 35;

%% Norm
mSignal  = enframe( vSignalNorm, hamming( iFrame ), iSolape );

vFrame = mSignal(idx,:);
vFrame = normalize( vFrame, 'zscore' );

% RPDE
[H_norm, rpd] = rpde( vFrame', m, tau );
bar(1:length(rpd),rpd,30)

str = ['$RPDE=', num2str( H_norm, 2 ), '$'];
TextLocation( str,'Interpreter','latex', 'Location', 'Best' )

name = 'RPDE_Norm'; 
axisIn  = [0 length(rpd) 0 0.6];
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )

%% Path
mSignal  = enframe( vSignalPath, hamming( iFrame ), iSolape );

vFrame = mSignal(idx,:);
vFrame = normalize( vFrame, 'zscore' );

% RPDE
[H_norm, rpd] = rpde( vFrame', m, tau );
bar(1:length(rpd),rpd,30)

str = ['$RPDE=', num2str( H_norm, 2 ), '$'];
TextLocation( str,'Interpreter','latex', 'Location', 'Best' )

name = 'RPDE_Path'; 
axisIn  = [0 length(rpd) 0 0.6];
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )