 function H = MelFilt(Nbank,NFFT,lopass,hipass,fs)
%function H = MelFilt(Nbank,NFFT,lopass,hipass,fs)
%
% Creates "Nbanks" mel-frequency triangular filters, H(Nbank,NFFT/2+1)
% First row of H corresponds to the lower CF filter.
%
% Nbank		- number of triangular filters in ]0..pi[
% NFFT		- number of FFT points [0..2*pi[
%               NOTE: H has points in [0..pi] : NFFT/2+1 points. 
% lopass	- lower cut-off of the first filter
% hipass	- upper cut-off of the last filter
% fs		- sampling frequency
%
% H may be directly applied to a FFT result:
% Y=fft(X);		  % with X=X(NFFT,Nframes): 1 column per frame
% Y=abs(Y(1:NFFT/2+1,:)); % abs in [0..pi]
% Z=H*Y;	% amplitude Filterbank result: Z=Z(Nbank,Nframes)
%		% (1 column (parameter vector) per frame)
%
% NOTE: according to HTK_V2.1, i.e. the filters are triangular on a mel scale,
%       not in the linear frequency scale.
 
% F. Perdigao
% Coimbra, Feb. 97

f=linspace(0,fs/2,NFFT/2+1);	%-- k=[0..NFFT/2] or [0..pi] rad/s

mlo=1127*log(1+lopass/700);
mhi=1127*log(1+hipass/700);
melCF=linspace(mlo,mhi,Nbank+2);	%-- corner points and CFs in mels
%fCF=700*(exp(mel/1127)-1);		%-- melCF in linear frequency
mel=1127*log(1+f/700);


%---- set the weigths from [0..pi] on the MEL scale -----------------

for k=1:Nbank,
 H(k,:)=triangf(mel,melCF(k+1),melCF(k),melCF(k+2));	%-- triangf(t,Tcenter,Tlow,Tupp)
end

%---- set the weigths from [0..pi] on the FREQUENCY scale -----------------
%for k=1:Nbank,
% H(k,:)=triangf(f,fCF(k+1),fCF(k),fCF(k+2));	%-- triangf(t,Tcenter,Tlow,Tupp)
%end



%----- to apply to an interval [0..2*pi[ :
%H=[H , fliplr(H(:,2:NFFT/2))];




