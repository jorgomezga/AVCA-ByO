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

idx = 13;

%% Calculation of parameters
[ ~, eParIntermedios ] = ...
    EM( vSignalPath, iFs, sModoPost, iNumParam, iFrame, iSolape, eOptions, iVerbosity );

eParIntermedios.sNombre =  '';
ePar.ePar = eParIntermedios;
[mParametros] = EM_ePar( ePar, vSignalPath, sModoPost, iNumParam, eOptions, iVerbosity);

bTuned = false;
[mParamHist] = histoParamEM_mod( eParIntermedios.vmMS, bTuned,...
    eParIntermedios.vMf, eParIntermedios.vAf, idx );

%% 
titulo  = '';
lineW   = 3;
fontS   = 35;
axisIn  = [-100 0 0 15e4];
bColor  = 1;
xLab = 'Energy (dB)';
yLab = 'Occurances';

name = 'CIL_Path'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS, bColor )