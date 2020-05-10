clear variables
close all
clc

addpath( genpath( '../../../../src' ) )
addpath( genpath( '../../../../External Toolboxes' ) )

addpath(  '../../' )

c2tdata = importdata( 'c2tdat_Path.mat' );

%% Si se calcula usando NL_TISEAN_d2 directamente
valNLTISEAN_d2 = importdata( 'HCTSAPath.mat' );
val_D2_mean = valNLTISEAN_d2.takens05_mean;

%% Si se calcula usando TS_COMPUTE
% valNLTISEAN_d2 = importdata( 'HCTSAPath.mat' );
% ID = strcmp( valNLTISEAN_d2.Operations.Name, 'NL_TISEAN_d2_1_10_0_takens05_mean' );
% val_D2_mean = valNLTISEAN_d2.TS_DataMat(ID);

xval = log10( c2tdata{1}(:,1) );
figure
for i=1:length(c2tdata)
    
    plot( log10( c2tdata{i}(:,1) ), c2tdata{i}(:,2), ...
    'LineWidth', 4 )
    hold on
end

% plot linea
plot( [xval(1), xval(end)], [val_D2_mean, val_D2_mean], '--k', ...
    'LineWidth', 4 );
str = ['$D_2$=',num2str(val_D2_mean)];
% TextLocation( str, 'Location', 'Best' )
valFin = xval(end)-0.1;
text(xval(1), valFin, str,...
    'Interpreter', 'LaTex', 'Fontsize', 35)

xLab    = '$\log{(\epsilon)}$';
yLab    = '$\mathcal{J}(\epsilon)$';
titulo  = '';
axisIn  = [-2.3 0 0 8];
lineW   = 3;
fontS   = 35;

name = 'D2_Path'; 
Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS )