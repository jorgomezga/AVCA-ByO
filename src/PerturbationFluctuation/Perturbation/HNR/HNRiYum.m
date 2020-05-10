function rHNRi=HNRiYum( Ci )

% Calculates the value of the HNR ratio for a speech segment corresponding to
% various pitch periods, collected in Ci, according to the Yumoto method
%
%   Yumoto, E., Gould, W. J., & Baer, T. (1982). Harmonics‐to‐noise ratio 
%   as an index of the degree of hoarseness. The journal of the Acoustical 
%   Society of America, 71(6), 1544-1550.
%
%   Yumoto, E., Sasaki, Y., & Okamura, H. (1984). Harmonics-to-noise ratio 
%   and psychophysical measurement of the degree of hoarseness. Journal of 
%   Speech, Language, and Hearing Research, 27(1), 2-6.
%
% Input parameters:
%   Ci:  matrix containing rows corresponding to pitch periods
%        of a segment of the voice signal, all of them with the duration
%        of the largest period (filled with zeros if necessary)
%
% Output parameters:
% rHNRi: the value of the HNR

rHNRi=0;
   
[n_periodos, ~] = size( Ci );

% An estimate of the harmonic component is obtained by taking the average of all
% the periods, since it makes possible to cancel the noise component (if the number
% of periods is sufficient).
Ca=mean(Ci);
% NOTA: la función mean calcula la FILA media de la matriz de entrada. Ca es un vector
%       fila de periodoMax elementos

H=n_periodos*sum(Ca.^2);

% The noise component for each pitch period can be calculated as the difference
% between said period and the average period.
N=0;
for i=1:n_periodos
   CiNoCeros=QuitaExtremos(Ci(i,:));
   L=length(CiNoCeros);
   Cr=CiNoCeros'-Ca(1:L);
   N=N+sum(Cr.^2);
end

if N~=0
   rHNRi=10*log10(H/N);
end