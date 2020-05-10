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


xLab    = '$\log{L}$';
yLab    = '${B}(L)$';
titulo  = '';
lineW   = 3;
fontS   = 35;
axisIn  = [0.7 3 -3 2];

%% Norm
mSignal  = enframe( vSignalNorm, hamming( iFrame ), iSolape );
vFrame = mSignal(idx,:);

[xpts, ypts] = fastdfa_core( vFrame' );

% Sort the intervals, and produce a log-log straight line fit
coeffs    = polyfit( log10(xpts), log10(ypts), 1);
vals      = polyval( coeffs, log10(xpts) );
alpha     = coeffs(1);
plot( log10(xpts), log10(ypts), 'r', log10(xpts), vals, 'k--', ...
    'LineWidth', 6 );

str = ['DFA=', num2str(alpha,2)];
TextLocation( str,'Interpreter','latex', 'Location', 'Best' )

name = 'DFA_Norm'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )

%% Path
mSignal  = enframe( vSignalPath, hamming( iFrame ), iSolape );
vFrame = mSignal(idx,:);

[xpts, ypts] = fastdfa_core( vFrame' );

% Sort the intervals, and produce a log-log straight line fit
coeffs    = polyfit( log10(xpts), log10(ypts), 1);
vals      = polyval( coeffs, log10(xpts) );
alpha     = coeffs(1);
plot( log10(xpts), log10(ypts), 'r', log10(xpts), vals, 'k--', ...
    'LineWidth', 6 );

str = ['DFA=', num2str(alpha,2)];
TextLocation( str,'Interpreter','latex', 'Location', 'Best' )

name = 'DFA_Path'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )