clear variables
clc
close all

addpath( genpath( '../../../../src' ) )
addpath( genpath( '../../../../External Toolboxes' ) )

sDir = '../../../../Audios';

% vSignalPath  = audioread( fullfile( sDir, '1424-a_n.wav' ) );
[vSignalPath, iFs]  = audioread( fullfile( sDir, 'cgra.wav' ) );

vSignalPath = normalize( vSignalPath, 'zscore' );

if exist( 'HCTSAPath.mat', 'file')
    delete( 'HCTSAPath.mat' )
end

idx       = 3;
iFrame    = ceil( 40e-3*iFs ); 
rSolape   = 0.5;
iSolape   = floor( (1 - rSolape)*iFrame );
sTipo     = 'M';
eOptions  = [];
iVerbosity = 1;

mSignal  = enframe( vSignalPath, hamming( iFrame ), iSolape );

%% Path
vFrame      = normalize( mSignal(idx,:), 'zscore' );

%lab = cellfun( @(x)['Op',num2str(x)], num2cell( [1:length(vFrame)] ), 'uni', 0 );
labels = {'frame'}; % data labels for each time series
keywords = {'signal1'}; % comma-delimited keywords for each time series
timeSeriesData = {vFrame}; % (a cell of vectors)

% Save these variables out to INP_test.mat:
save('INP_testPath.mat','timeSeriesData','labels','keywords');

% Initialize a new hctsa analysis using these data and the default feature library:
TS_init('INP_testPath.mat', 'INP_mops.txt', 'INP_ops.txt', false, 'HCTSAPath.mat');

%% Option 1: Focusing on what we need
out = NL_TISEAN_d2( vFrame', 1, 10, 0, 'c2tdat_Path.mat' );
% out = NL_TISEAN_d2( vFrame', 14, 10, 20 )
save( 'HCTSAPath.mat', 'out' )

%% Option 2: Using the whole toolbox
% NL_TISEAN_D2 in the hctsa tolbox has been modified, adding the prefix
% dbstop in NL_TISEAN_d2 at 205, i.e:
% dbstop c2tdat = SUB_readTISEANout(s,maxm,'#m=',2);
% This value is manually saved and this is the file to be plotted
% TS_compute(false,[],[],'missing','HCTSAPath.mat',true)