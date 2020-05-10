clear variables
close all

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

iFrame    = ceil( 150e-3*iFs );
rSolape   = 0;
iSolape   = floor( (1 - rSolape)*iFrame );

bSmooth   = 1;
sNormaliza     = 'line';

xLab    = 'Time (s)';
yLab    = 'Jitter (%)';
titulo  = '';
axisIn  = [0 2 0 30];
lineW   = 3;
fontS   = 30;

iNumPoints = 150;

%% Norm
% For the calculation of the disturbances and tremor some information is needed that
% is obtained from the pitch analysis
[vPPS, ~] = ParamPitch( vSignalNorm, iFs, 1, length( vSignalNorm ), 1 );

% jitter envelope
xVal = linspace( 0, length(vSignalNorm)/iFs, iNumPoints  );
[~, vJitter]=jitter_env( vPPS, iNumPoints );
stem( xVal, vJitter )

name = 'Jitter_Norm';
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )

%% Path
% For the calculation of the disturbances and tremor some information is needed that
% is obtained from the pitch analysis
[vPPS, ~] = ParamPitch( vSignalPath, iFs, 1, length( vSignalPath ), 1 );

% jitter envelope
xVal = linspace( 0, length(vSignalPath)/iFs, iNumPoints  );
[~, vJitter]=jitter_env( vPPS, iNumPoints );
stem( xVal, vJitter )

name = 'Jitter_Path';
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )