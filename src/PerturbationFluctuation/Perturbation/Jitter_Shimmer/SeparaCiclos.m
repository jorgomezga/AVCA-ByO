function vCiclos=SeparaCiclos( vSignal, iInicio, iFinal, i_P ) 

% Separates the cycles of a known voice signal given the moments in which
%    each cycle begins.
%
% Input parameters:
%   vSignal:    Speech signal
%   iInicio:    is the first sample of the section to be analyzed
%   iFinal:     is the last sample of the section to be analyzed
%   i_P:        the start positions of the different vSignal cycles.
%               (This parameter is obtained from applying the PicosPos or PicosNeg function
%               (Boyanov), or to use PicosPosKas and PicosNegKas from Kasuya).
% Output parameters
%   vCycles:    is a column array whose rows are the voice segments corresponding to consecutive cycles.
%               All segments will have the longest cycle length, filling in with zeros if necessary.

iNumPicos=length( i_P );

% Only cycles between two consecutive peaks will be considered,
% discarding the initial voice segments (up to the first peak) and iFinal (from
% the last peak);
iNumCiclos=iNumPicos-1;

% Duration in samples of the largest cycle
iDurCilcoMax=0;
vDurCiclo= zeros( 1, iNumCiclos); 
for i=1:iNumCiclos
    vDurCiclo(i)=i_P(i+1)-i_P(i);
    if vDurCiclo(i) > iDurCilcoMax
        iDurCilcoMax = vDurCiclo( i );
    end
end

voz=vSignal';

for i=1:iNumCiclos
    vCiclos(i,:)=[voz( i_P(i):i_P(i+1)-1 ) zeros( 1 ,iDurCilcoMax-vDurCiclo(i) )];
end

if nargout == 0
    % Representamos los 5 primeros ciclos
    for i=1:4, subplot(4,1,i); plot(vCiclos(i,:)); end
    title ('componente armónica');
end