clear variables
close all
clc

addpath(genpath('../../'))
addpath(genpath('../../../libs'))

sDir = '../../Audios';
[vSignalNorm, iFs]  = audioread( fullfile( sDir, 'asra.wav' ) );
vSignalPath  = audioread( fullfile( sDir, 'cgra.wav' ) );

%% Parameters
sTipo='MCYFOPA'; % All features are calculated, using Hamming windowing
iSubMod=26;
iFrame=pow2(floor(log2(0.150*iFs))); % Default: 150 ms
iShift=floor(iFrame/2);
eOptions=[];
iVerbosity=0;

%% Norm
mFeaturesNorm = MS_features( vSignalNorm, iFs, sTipo,...
    iSubMod, iFrame, iShift, eOptions, iVerbosity);

%% Path
mFeaturesPath = MS_features( vSignalPath, iFs, sTipo,...
    iSubMod, iFrame, iShift, eOptions, iVerbosity);
