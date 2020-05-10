clear variables
close all
clc

addpath( genpath( '../src' ) )
addpath( genpath( '../External Toolboxes' ) )

sDir = '../Audios';

[vSignalNorm, iFs]  = audioread( fullfile( sDir, 'asra.wav' ) );
vSignalNorm = normalize( vSignalNorm, 'zscore' );

vSignalPath  = audioread( fullfile( sDir, 'cgra.wav' ) );
vSignalPath = normalize( vSignalPath, 'zscore' );

%% Parameters
iNumberOfFrames = 100;

%% Norm
sType = 'Boyanov';
[~, vFoNorm_Boyanov] = CalculatePitch( vSignalNorm, iFs, iNumberOfFrames, sType );
sType = 'Rabiner_StaticClip';
[~, vFoNorm_Rabiner_StaticClip] = CalculatePitch( vSignalNorm, iFs, iNumberOfFrames, sType );
sType = 'Rabiner_LPC';
[~, vFoNorm_Rabiner_LPC] = CalculatePitch( vSignalNorm, iFs, iNumberOfFrames, sType );
sType = 'Rabiner_Cepstral';
[~, vFoNorm_Rabiner_Cepstral] = CalculatePitch( vSignalNorm, iFs, iNumberOfFrames, sType );
sType = 'Kasuya';
[~, vFoNorm_Kasuya] = CalculatePitch( vSignalNorm, iFs, iNumberOfFrames, sType );


figure(1)

xTime = linspace( 0, length( vSignalNorm )/iFs, iNumberOfFrames );
plot( xTime, vFoNorm_Boyanov, 'r' ), hold on
plot( xTime, vFoNorm_Rabiner_StaticClip, 'y' ), hold on
plot( xTime, vFoNorm_Rabiner_LPC, 'r' ), hold on
plot( xTime, vFoNorm_Rabiner_Cepstral, 'm--' ), hold on


xTime = linspace( 0, length( vSignalNorm )/iFs, length( vFoNorm_Kasuya ) );
plot( xTime, vFoNorm_Kasuya, 'g' ), hold on
legend( {'Boyanov', 'Rabiner StaticClip', 'Rabiner LPC', 'Rabiner Cepstral', 'Kasuya' } )

%% Path
sType = 'Boyanov';
[~, vFoPath_Boyanov] = CalculatePitch( vSignalPath, iFs, iNumberOfFrames, sType );
sType = 'Rabiner_StaticClip';
[~, vFoPath_Rabiner_StaticClip] = CalculatePitch( vSignalPath, iFs, iNumberOfFrames, sType );
sType = 'Rabiner_LPC';
[~, vFoPath_Rabiner_LPC] = CalculatePitch( vSignalPath, iFs, iNumberOfFrames, sType );
sType = 'Rabiner_Cepstral';
[~, vFoPath_Rabiner_Cepstral] = CalculatePitch( vSignalPath, iFs, iNumberOfFrames, sType );
sType = 'Kasuya';
[~, vFoPath_Kasuya] = CalculatePitch( vSignalPath, iFs, iNumberOfFrames, sType );

figure(2)

xTime = linspace( 0, length( vSignalPath )/iFs, iNumberOfFrames );
plot( xTime, vFoPath_Boyanov, 'r' ), hold on
plot( xTime, vFoPath_Rabiner_StaticClip, 'y' ), hold on
plot( xTime, vFoPath_Rabiner_LPC, 'r' ), hold on
plot( xTime, vFoPath_Rabiner_Cepstral, 'm--' ), hold on

xTime = linspace( 0, length( vSignalPath )/iFs, length( vFoPath_Kasuya ) );
plot( xTime, vFoPath_Kasuya, 'g' ), hold on

legend( {'Boyanov', 'Rabiner StaticClip', 'Rabiner LPC', 'Rabiner Cepstral', 'Kasuya' } )