clear variables
close all
clc

addpath( genpath( '../../../../src' ) )
addpath( genpath( '../../../../External Toolboxes' ) )

sDir = '../../../../Audios';

% vSignalNorm  = audioread( fullfile( sDir, '1-a_n.wav' ) );
[vSignalNorm, iFs]  = audioread( fullfile( sDir, 'asra.wav' ) );
vSignalNorm = normalize( vSignalNorm, 'zscore' );

idx       = 3;
iFrame    = ceil( 40e-3*iFs ); 
rSolape   = 0.5;
iSolape   = floor( (1 - rSolape)*iFrame );


mSignal  = enframe( vSignalNorm, hamming( iFrame ), iSolape );
vFrame    = mSignal(idx,:);

iTao = 5;
iDim = 4;
mAtractor   = embeb( vFrame, iDim, iTao );  

%% Norm
iExclude    = ceil( 0.1*length( vFrame ) );  % in case the query points are taken out of the pointset, 
                                             % exclude specifies a range of indices which are omitted fromsearch.
% Largest Lyapunov
iMaxdelta_n = 30;  % Maximo de iteraciones, para el calculo de delta_n
iNumVecinos = 5;  % Number of nearest neighbors to compute

%x = largelyap(pointset, query_indices, taumax, maximal_neighbours, k_exclude)
l    = largelyap_R( mAtractor, 1:length( mAtractor )-iMaxdelta_n, iMaxdelta_n, iNumVecinos, iExclude );
LLE  = polyfit( 1:iMaxdelta_n, l(1:iMaxdelta_n)', 1 );
iLLE = LLE(1);

xpts    = 1:iMaxdelta_n;
vals    = polyval( LLE, xpts );
plot( xpts, l(1:iMaxdelta_n)', 'r', xpts, vals, 'k--', ...
    'LineWidth', 6 )

str = ['LLE=', num2str(iLLE,1)];
TextLocation( str, 'Location', 'best', 'FontSize', 14, 'Interpreter', 'LaTex' ); 

xLab    = '$\Delta{t}$';
yLab    = '$Q(\Delta{t})$';
titulo  = '';
axisIn  = [1 30 -3 0];
lineW   = 3;
fontS   = 35;

name = 'LLE_Norm'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )
