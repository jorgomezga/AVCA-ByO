clear variables
clc
close all

addpath( genpath( '/mnt/A0D63A36D63A0D52/QoV-ByO' ) )


% addpath( genpath( '/home/jorge/Documentos/Repositorio' ) )
% addpath( '/mnt/data/Dropbox/ResultadosThesis/AuxFunctions/axlabel' )

% sIn = '/mnt/data/BasesDatos/BD_HUPA/BDHUPAseg_UoC_16_a_69_rev2/PorCondicion/Normal/aaa.wav';
% sIn = '/mnt/data/BasesDatos/BD_HUPA/BDHUPAseg_UoC_16_a_69_rev2/PorCondicion/Pathol/araapath.wav';


% sIn = '/mnt/data/BasesDatos/KayElemetrics/PorCondicion/normal/AXH1NAL.WAV';
% sIn = '/mnt/data/BasesDatos/KayElemetrics/PorCondicion/pathol/DVD19AN.WAV';
% sIn = '/mnt/data/BasesDatos/KayElemetrics/PorCondicion/pathol/ALB18AN.WAV';

sIn = '/mnt/A0D63A36D63A0D52/Databases/Kay/Ah/Normal/AXH1NAL.WAV';
% sIn = '/mnt/A0D63A36D63A0D52/Databases/Kay/Ah/Pathol/DVD19AN.WAV';
[vSignal, iFs] = audioread( sIn );
    
iFrame   = ceil( 40e-3*iFs );
iSolape  = ceil( 0.5*iFrame );
sTipo    = 'R';
iNumCoef = 14;
eOptions = [];
iVerbosity = 1;

vSignal = vSignal(1:25000);
MFCCs( vSignal, iFs, sTipo, iNumCoef, iFrame,...
    iSolape, eOptions, iVerbosity )
    axis([-0.2 0.2])