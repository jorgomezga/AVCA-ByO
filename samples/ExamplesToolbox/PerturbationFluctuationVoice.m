clear variables
close all
clc

addpath(genpath('../../'))
addpath(genpath('../../../libs'))

sDir = '../../Audios';
[vSignalNorm, iFs]  = audioread( fullfile( sDir, 'asra.wav' ) );
vSignalPath  = audioread( fullfile( sDir, 'cgra.wav' ) );

iFrame      = ceil( 40e-3*iFs ); 
iOverlap    = floor( 0.5*iFrame );

%% Norm
[vJitterNorm, vShimmerNorm]  = JitterShimmer( vSignalNorm, iFs );
vFluctuationNorm = Fluctuation( vSignalNorm, iFs );
vAdditiveNoiseNorm = AdditiveNoise( vSignalNorm, iFs, iFrame, iOverlap );

%% Path
[vJitterPath, vShimmerPath]  = JitterShimmer( vSignalPath, iFs );
vFluctuationPath = Fluctuation( vSignalPath, iFs );
vAdditiveNoisePath = AdditiveNoise( vSignalPath, iFs, iFrame, iOverlap );