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

iFrame    = ceil( 100e-3*iFs ); 
rSolape   = 0.8;
iSolape   = floor( (1 - rSolape)*iFrame );

bSmooth   = 1;
sNormaliza     = 'line';

xLab    = 'Time (s)'; 
yLab    = 'Shimmer (\%)'; 
titulo  = '';
axisIn  = [0 1.2 0 80];
lineW   = 3;
fontS   = 35;

iNumPoints = 150;

%% Norm
% For the calculation of the disturbances and tremor some information is needed that
% is obtained from the pitch analysis
[~, vPAS] = ParamPitch( vSignalNorm, iFs, 1, length( vSignalNorm ), 1 );

% jitter envelope
xVal = linspace( 0, length(vSignalNorm)/iFs, iNumPoints  );
[~, rShim] = shimmer_env( vPAS, iNumPoints );
stem( xVal, rShim )

name = 'Shimmer_Norm'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )
    
%% Path
% For the calculation of the disturbances and tremor some information is needed that
% is obtained from the pitch analysis
[~, vPAS] = ParamPitch( vSignalPath, iFs, 1, length( vSignalPath ), 1 );

% jitter envelope
xVal = linspace( 0, length(vSignalPath)/iFs, iNumPoints  );
[~, rShim] = shimmer_env( vPAS, iNumPoints );
stem( xVal, rShim )

name = 'Shimmer_Path'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )