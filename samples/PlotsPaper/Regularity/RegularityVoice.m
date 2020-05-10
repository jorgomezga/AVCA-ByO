clear variables
close all
clc

addpath( genpath( '/mnt/A0D63A36D63A0D52/QoV-ByO/' ) )


sDirNorm = '/mnt/A0D63A36D63A0D52/Databases/Saarbruecken/Norm';
sDirPath = '/mnt/A0D63A36D63A0D52/Databases/Saarbruecken/Path';

iFs       = 50e3;
iFrame    = ceil( 40e-3*iFs ); 
rSolape   = 0.5;
iSolape   = floor( (1 - rSolape)*iFrame );
sTipo     = 'M';
eOptions  = [];
iVerbosity = 1;



% vSignal = audioread( '/mnt/A0D63A36D63A0D52/Databases/Kay/Ah/Normal/AXH1NAL.WAV' );
% vSignal = audioread( '/mnt/A0D63A36D63A0D52/Databases/Kay/Ah/Pathol/DVD19AN.WAV' );

%%%% NO
% vSignal  = audioread( fullfile( sDirPath, '1424-a_n.wav' ) );
% vSignal  = audioread( fullfile( sDirNorm, '1-a_n.wav' ) );
%%%%

%% Norm
vSignal  = audioread( fullfile( sDirNorm, '1-a_n.wav' ) );
mSignal  = Enventanado( vSignal, iFs, iFrame, iSolape, sTipo, eOptions, iVerbosity );

dim = 8;
tau = 20;
iAlpha = 0.35;

rApEn = zeros( 1, size( mSignal, 1 ) );
rSampEn = zeros( 1, size( mSignal, 1 ) );
rFuzzyEn = zeros( 1, size( mSignal, 1 ) );
rGSampEn = zeros( 1, size( mSignal, 1 ) );
rmSampEn = zeros( 1, size( mSignal, 1 ) );

xNorm = 1:size( mSignal, 1 );

for j=1:size( mSignal, 1 )

    vFrame = mSignal(j,:);
    vFrame = Normaliza( vFrame, 'unomenosuno' );    
    rParam = std( vFrame )*iAlpha;
    
    % Reconstruccion
    % dim = m + 1
    mAtractorEntropyp1 = embeb( vFrame, dim+1, tau );
    % dim = m
    mAtractorEntropy = embeb( vFrame, dim, tau );
    
    [rApEn(j),rSampEn(j),rFuzzyEn(j),rGSampEn(j),rmSampEn(j)] = ...
        CalculateEntropies( mAtractorEntropy, mAtractorEntropyp1, rParam, 1 ); 
    
end

%% Path
vSignal  = audioread( fullfile( sDirPath, '1424-a_n.wav' ) );
mSignal  = Enventanado( vSignal, iFs, iFrame, iSolape, sTipo, eOptions, iVerbosity );

dim = 8;
tau = 20;

rApEnPath = zeros( 1, size( mSignal, 1 ) );
rSampEnPath = zeros( 1, size( mSignal, 1 ) );
rFuzzyEnPath = zeros( 1, size( mSignal, 1 ) );
rGSampEnPath = zeros( 1, size( mSignal, 1 ) );
rmSampEnPath = zeros( 1, size( mSignal, 1 ) );
for j=1:size( mSignal, 1 )

    vFrame = mSignal(j,:);
    vFrame = Normaliza( vFrame, 'unomenosuno' );    
    rParam = std( vFrame )*iAlpha;
    
    % Reconstruccion
    % dim = m + 1
    mAtractorEntropyp1 = embeb( vFrame, dim+1, tau );
    % dim = m
    mAtractorEntropy = embeb( vFrame, dim, tau );
    
    [rApEnPath(j),rSampEnPath(j),rFuzzyEnPath(j),rGSampEnPath(j),rmSampEnPath(j)] = ...
        CalculateEntropies( mAtractorEntropy, mAtractorEntropyp1, rParam, 1 ); 
end

figure
plot( xNorm, rApEn, xNorm, rSampEn, xNorm, rFuzzyEn, xNorm, rGSampEn, xNorm, rmSampEn,...
    'LineWidth', 6 )

axis( [1, length( xNorm ), 0, .6] )
ylabel( 'Regularity', 'Interpreter','LaTex' )
xlabel( 'window','Interpreter','LaTex' )

legend( {'ApEn', 'SampEn', 'FuzzyEn', 'GSampEn', 'mSampEn'}, 'Fontsize', 40 )

set(gca,'fontsize',35)
pretty_xyplot

%%%%%
figure
xPath = 1:size( mSignal, 1 );
plot( xPath, rApEnPath, xPath, rSampEnPath, xPath, rFuzzyEnPath, xPath,...
    rGSampEnPath, xPath, rmSampEnPath, 'LineWidth', 6 )

axis( [1, length( xNorm ), 0, .6] )
ylabel( 'Regularity', 'Interpreter','LaTex' )
xlabel( 'window','Interpreter','LaTex' )

legend( {'ApEn', 'SampEn', 'FuzzyEn', 'GSampEn', 'mSampEn'}, 'Fontsize', 40 )

set(gca,'fontsize',35)
pretty_xyplot