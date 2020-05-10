function [MelFBank] = Create_MelFrequencyFilterBank(fs, Nfft, Nfilt)
%Create_MelFrequencyFilterBank 
%
%[MelFBank] = Create_MelFrequencyFilterBank(fs, Ndft, Nfilt)
%creates a Mel-scale filter bank to be used with a signal, whose sampling frequency 
%is 'fs' and which is transformed with DFT having 'Ndft' frequency bins. 
%The function creates 'Nfilt' number of triangular filters,
%which are created to be equally spaced in Mel-frequency scale,
%and returns them as rows of matrix 'MelFBank' whose size is 'Nfilt' x 'Nfft/2'.
%
%Create_MelFrequencyFilterBank(fs, Ndft, Nfilt) 
%Without output arguments, the function plots the created Mel-scale filter bank.

maxf = fs/2;                                    % The maximal frequency used
maxmelf = 2595*log10(1+maxf/700);               % The maximal Mel-frequency -value  
edgemelfs = (0:(Nfilt+1))/(Nfilt+1) * maxmelf;  % All triangle-filter edge Mel-frequency -values
edgefrqs = 700*(10.^(edgemelfs/2595)-1);        % All triangle-filter edge frequencies ( in normal frequency scale )
edgeDFTbins = round(edgefrqs/maxf*(Nfft/2));    % All triangle-filter edge DFTbins customized to Nfft
edgeDFTbins(1) = 1;                             % Just correcting the first bin not to be 0.

if edgeDFTbins(2) == 0,
    error('Can not create so many filters for this Ndft!')
end

MelFBank = zeros(Nfilt,Nfft/2);                 % Reserving memory for the filterbank. Each row of 'MelFBank'
                                                % corresponds to one triangle-filter.
for n=1:Nfilt,
    l = edgeDFTbins(n);                         % start index for the lower edge
    c = edgeDFTbins(n+1);                       % center index for the triangle
    h = edgeDFTbins(n+2);                       % end index for the higher edge

    NbinsUpSlope = c-l;                         % Number of DFT-points in lower edge
    NbinsDownSlope = h-c;                       % Number of DFT-ponts in higher edge

    MelFBank(n,l:c) = (0:NbinsUpSlope)/NbinsUpSlope;        % Create the lower (frequency) slope  (going up) ...
    MelFBank(n,c:h) = (NbinsDownSlope:-1:0)/NbinsDownSlope; % Create the higher (frequency) slope (going down) ...
end                                                         %  .. of filter number n.

%_____
if ~nargout
figure('Name','Mel Frequency Filter Bank')
problem4fig = gcf;
set(problem4fig,'position',[1 700 560 420]);
plot(MelFBank')
axis([0, Nfft/2,  0,  1.1])
title('Mel-filterbank')
end