clear variables
close all
clc

addpath( genpath( '/home/jorge/qov-byo' ) )

sDir = '/home/jorge/Dropbox/Toolbox_Part3/Audios';
vSignalNorm  = audioread( fullfile( sDir, '1-a_n.wav' ) );
vSignalPath  = audioread( fullfile( sDir, '1424-a_n.wav' ) );

%% Parameters
iFs       = 50e3;
iFrame    = ceil( 40e-3*iFs ); 
rSolape   = 0.5;
iSolape   = floor( (1 - rSolape)*iFrame );

%% Norm


% [vParametros, caNombres] = PerturbationFluctuation( vSignalNorm, iFs );
[vParametros, caNombres] = PerturbationFluctuation( vSignalPath, iFs );

% Use the voicebox toolbox
mSignal = enframe( vSignalNorm, hamming( iFrame ), iSolape );

for j=1:size( mSignal, 1 )

    vFrame = mSignal(j,:);
    vFrame = normalize( vFrame, 'zscore' );    
    
end
