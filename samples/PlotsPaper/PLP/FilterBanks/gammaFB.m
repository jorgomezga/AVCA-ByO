  function [H,BGAMMA,AGAMMA,CF]=gammaFB(W)
% function [H,BGAMMA,AGAMMA,CF]=gammaFB(W);
%
% Definition of a gammatone filterbank model for 
% 8kHz sampling frequency. 35 channels. 
% Uses function MakeERBFilters from Slaney's Auditory Toolbox


fs=8000;
[B,A,CF]=MakeERBFilters(fs,35,200);

%f=logspace(log10(50),log10(4000),500);
%W=2*pi*f/fs;

for k=1:35,
 H(k,:)=freqz(B(k,:),A(k,:),W);
end
%semilogx(f,db(H))


BGAMMA=B; AGAMMA=A;
%save GAMA BGAMMA AGAMMA

