clear variables
close all
clc

addpath(genpath('../../'))
addpath(genpath('../../../libs'))

sDir = '../../Audios';
[vSignalNorm, iFs]  = audioread( fullfile( sDir, 'asra.wav' ) );
vSignalPath  = audioread( fullfile( sDir, 'cgra.wav' ) );

%% Parameters
iFrame    = ceil( 40e-3*iFs ); 
iSolape   = floor( 0.5*iFrame );
sNorm     = 'chebychev';

dim       = 8;
tau       = 20;
iAlpha    = 0.35;

%% Norm
% Use the voicebox toolbox
mSignal = enframe( vSignalNorm, hamming( iFrame ), iSolape );

rApEn = zeros( 1, size( mSignal, 1 ) );
rSampEn = zeros( 1, size( mSignal, 1 ) );
rFuzzyEn = zeros( 1, size( mSignal, 1 ) );
rGSampEn = zeros( 1, size( mSignal, 1 ) );
rmSampEn = zeros( 1, size( mSignal, 1 ) );

xNorm = 1:size( mSignal, 1 );
for j=1:size( mSignal, 1 )

    vFrame = mSignal(j,:);
    vFrame = normalize( vFrame, 'zscore' );    
    rParam = std( vFrame )*iAlpha;
    
    % Reconstruction
    % dim = m + 1
    mAtractorEntropyp1 = embeb( vFrame, dim+1, tau );
    % dim = m
    mAtractorEntropy = embeb( vFrame, dim, tau );
        
    [rApEn(j),rSampEn(j),rmSampEn(j),rGSampEn(j),rFuzzyEn(j)] =...
        CalculateRegularity( mAtractorEntropyp1, mAtractorEntropy, rParam, sNorm, true );
end

figure
plot( xNorm, rApEn, xNorm, rSampEn, xNorm, rFuzzyEn, xNorm, rGSampEn, xNorm, rmSampEn,...
    'LineWidth', 6 )
title('Normophonic')

axis( [1, length( xNorm ), 0, .35] )
ylabel( 'Regularity', 'Interpreter','LaTex' )
xlabel( 'window','Interpreter','LaTex' )
legend( {'ApEn', 'SampEn', 'FuzzyEn', 'GSampEn', 'mSampEn'}, 'Fontsize', 20 )
set(gca,'fontsize',20)

%% Path
% Use the voicebox toolbox
mSignal = enframe( vSignalPath, hamming( iFrame ), iSolape );

rApEnPath = zeros( 1, size( mSignal, 1 ) );
rSampEnPath = zeros( 1, size( mSignal, 1 ) );
rFuzzyEnPath = zeros( 1, size( mSignal, 1 ) );
rGSampEnPath = zeros( 1, size( mSignal, 1 ) );
rmSampEnPath = zeros( 1, size( mSignal, 1 ) );

for j=1:size( mSignal, 1 )

    vFrame = mSignal(j,:);
    vFrame = normalize( vFrame, 'zscore' );   
    rParam = std( vFrame )*iAlpha;
    
    % Reconstruction
    % dim = m + 1
    mAtractorEntropyp1 = embeb( vFrame, dim+1, tau );
    % dim = m
    mAtractorEntropy = embeb( vFrame, dim, tau );
    
    [rApEnPath(j),rSampEnPath(j),rmSampEnPath(j),rGSampEnPath(j),rFuzzyEnPath(j)] =...
        CalculateRegularity( mAtractorEntropyp1, mAtractorEntropy, rParam, sNorm, true );
end

%%%%%
figure
xPath = 1:size( mSignal, 1 );
plot( xPath, rApEnPath, xPath, rSampEnPath, xPath, rFuzzyEnPath, xPath,...
    rGSampEnPath, xPath, rmSampEnPath, 'LineWidth', 6 )
title('Dysphonic')

axis( [1, length( xNorm ), 0, .35] )
ylabel( 'Regularity', 'Interpreter','LaTex' )
xlabel( 'window','Interpreter','LaTex' )
legend( {'ApEn', 'SampEn', 'FuzzyEn', 'GSampEn', 'mSampEn'}, 'Fontsize', 20 )
set(gca,'fontsize',20)