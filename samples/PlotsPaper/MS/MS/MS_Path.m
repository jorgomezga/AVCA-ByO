clear variables
close all
clc

addpath( '../../' )
addpath( genpath( '../../../../src' ) )
addpath( genpath( '../../../../External Toolboxes' ) )

sDir = '../../../../Audios';

% [vSignalPath, iFs]  = audioread( fullfile( sDir, '1424-a_n.wav' ) );
[vSignalPath, iFs]  = audioread( fullfile( sDir, 'cgra.wav' ) );

vSignalPath = normalize( vSignalPath, 'zscore' );


%Parametros Ventana
iTraslape = 0.5;
iFrame    = ceil( 180e-3*iFs );
iSolape   = ceil( iTraslape*iFrame );

%Procesado
sModoPrep = 'nHQ';
sModoPost = 'XCFPA';
%Centroids (normalized respecting the first centroid),
%Contrast & Homogeneity and Histogram parameters

%Parametros EMod
eOptions.EMiFmodMax=200;
iNumParam  = 12;
iVerbosity = 0;

%% Calculation of parameters
[ ~, eParIntermedios ] = ...
    EM( vSignalPath, iFs, sModoPost, iNumParam, iFrame, iSolape, eOptions, iVerbosity );

eParIntermedios.sNombre =  '';
ePar.ePar = eParIntermedios;
[mParametros] = EM_ePar( ePar, vSignalPath, sModoPost, iNumParam, eOptions, iVerbosity);


% for idx = 1:size( eParIntermedios.vmMS, 3 )

idx = 13;

temp = 20*log10( fftshift( abs( eParIntermedios.vmMS(:,:,idx) ), 2 ) );
figure
if 	~all( diff( eParIntermedios.vAf ) == min( diff( eParIntermedios.vAf ) ) )
    % This is the case where the subband center frequencies are
    % non-uniformly spaced. The following code interpolates P along the
    % acoustic-frequency axis using a base granularity set by STEPSIZE.
    % As a result the modulation spectra of broad subbands will appear
    % as thicker rows in the displayed joint-frequency plot.
    stepsize = min( 50, min( diff( eParIntermedios.vAf ) / 2 ) );
    afreqs2 = eParIntermedios.vAf(1):stepsize:eParIntermedios.vAf(end);
    temp = interp1( eParIntermedios.vAf, temp, afreqs2, 'nearest' );
    imagesc( eParIntermedios.vMf - eParIntermedios.data.modfs/2, afreqs2, temp );
else
    % Uniformly-spaced subbands are much easier to plot.
    imagesc( mfreqs - eParIntermedios.data.modfs/2, eParIntermedios.vAf, temp );
end

if strcmpi( eParIntermedios.data.demodparams{1}, 'harm' ) ||...
        strcmpi( eParIntermedios.data.demodparams{1}, 'harmcog' )
    axislabel = 'Harmonic number';
else
    axislabel = 'Acoustic frequency (Hz)';
end

%% 
xLab = 'Modulation frequency (Hz)';
yLab = axislabel;
axis xy

titulo  = '';
lineW   = 3;
fontS   = 35;
axisIn  = '';
bColor  = 1;

name = 'MS_Path'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS, bColor )
% end