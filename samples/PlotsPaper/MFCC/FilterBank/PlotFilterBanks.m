
addpath( genpath( '/mnt/A0D63A36D63A0D52/QoV-ByO' ) )


iLowFreq=0;
iHighFreq=0.5;
iFs = 16000;
iNFFT = 1024;


iNumFilters=floor(3*log(iFs));

melbankm(iNumFilters, iNFFT, iFs, iLowFreq, iHighFreq)


% n = 1024;
% p = 30;
% fs = 16000;
% x=filtbankm(p,n,fs,0,fs/2,'m');            % n is the fft length, p is the number of filters wanted
% 
% %  (c) Plot the calculated filterbanks
% 
% plot((0:floor(n/2))*fs/n,filtbankm(p,n,fs,0,fs/2,'m')')   % fs=sample frequency

%  (d) Plot the filterbanks

% filtbankm(p,n,fs,0,fs/2,'m');