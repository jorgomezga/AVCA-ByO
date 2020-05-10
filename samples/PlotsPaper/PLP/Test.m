
addpath( genpath( '/home/jorge/Documentos/Repositorio' ) ) 

addpath( genpath( '/mnt/A0D63A36D63A0D52/QoV-ByO' ) )


% iLowFreq=0;
% iHighFreq=0.5;
% iFs = 16000;
% iNFFT = 1024;
% 
% 
% iNumFilters=floor(3*log(iFs));
% 
% melbankm(iNumFilters, iNFFT, iFs, iLowFreq, iHighFreq, 'b')


N = 8000;
Nband = 19;
fs = 8000;

[H,E]=PLPfilt(N,Nband,fs);

H = full(H);

plot(H','b');
    set(gca,'xlim',[0 N]);
    xlabel(['Frequency (' xticksi 'Hz)']);
    