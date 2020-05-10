  function [c,P]=JRasta(x,N,M,Nband,fs,p,Q,Jah)
% function [c,P]=JRasta(x,N,M,Nband,fs,p,Q,Jah);
%
% LIN-LOG RASTA-PLP cepstral analysis of vector x
%
% Performs J-RASTA analysis on vector x, ...
% with frames of N samples, every M samples, with Nband filters.
%
% x: signal
% N: frame length
% M: frame rate in samples
% Nband: num. of (effective) critical bands
% fs: sampling frequency
% Jah: constant for non-linearity: y=log(1+Jah*x) (default: Jah=1e-6)
% c: vector with Q cepstral coeficients (including c0). 
% 
% NOTE: N,M,Nband,p,Q must be integers!
% NOTE: This function computs cepstrum according to 
%       Morgan's RASTA_2_2 C-code implementation.
%       [www.icsi.berkley.edu/ftp/global/pub/speech/rasta/]
% Ref.: "Rasta-PLP Speech Analysis Technique", H. Hermansky et al, ICASSP'92, pp. I-121.
% see also PLP.m PLPfilt.m

% Fernando Perdigao, Dec. 97
% fp@it.uc.pt


if nargin<7, error('too few parameters...'); end
if p>Nband, error('LP order must be less then Nband'); end
if nargin<8, Jah=1e-6; end

lift_exp=0.6;	%--- for liftering as h(n)=(n^lift_exp), n=1:Q
pole=0.94;	%--- for IIR filtering

%---------------------------------------------------------------------------

[X,nSamples,Nfft]=frames(x,N,M); %--- block into frames (with hamming window)

X=abs(fft(X));
X=X(1:Nfft/2+1,:).^2;	%-- range: [0..pi]

[H,Elp]=PLPfilt(Nfft/2+1,Nband,fs); %--- PLP filters and preemphasis factors

logP=log(1+ Jah*H*X);	 %--- Lin-Log power spectrum

%--- to behave exactly as in Morgan's RASTA_2_2 
%--- [www.icsi.berkley.edu/ftp/global/pub/speech/rasta] 

%rasta=filter((0.2:-0.1:-0.2),1,logP');	%--- fir filtering
%rasta(1:4,:)=zeros(4,Nband);		%--- not allow great change at begining
%rasta=filter(1,[1,-pole],rasta)';	%--- iir filtering

%---- for MatLab v4 -----
rasta=zeros(Nband,nSamples);
b=(0.2:-0.1:-0.2); a=[1,-pole];
for k=1:Nband,
 y=filter(b,1,logP(k,:));
 y(1:4)=zeros(1,4);
 rasta(k,:)=filter(1,a,y);
end


%--- another way (simpler): doesn't need to eliminate 
%--- the 4 initial frames ------
%rasta=deltas(logP)/sqrt(10);
%rasta=filter(1,[1,-pole],rasta')';





nlspec=exp(rasta);
nlspec=( nlspec.*Elp(:,ones(1,nSamples)) ).^0.33;


 %--- autocorrelation: ------
 R=real(ifft([nlspec(1,:);nlspec;nlspec(Nband,:);flipud(nlspec)]));

 %--- LP analysis ----------
 a=zeros(p+1,nSamples);
 MSE=zeros(1,nSamples);

  for k=1:nSamples,
    r=R(1:p+1,k);
    a(:,k) = [ 1; -toeplitz(r(1:p))\r(2:p+1)];	%-- LP coeficients
    MSE(k)=r'*a(:,k);				%--- MS error
  end

%----- cepstrum ----------
c=lpc2ceps(a,Q-1);

%--- liftering ---------
hlift=(1:Q-1)'.^lift_exp;
c=[log(MSE);c.*hlift(:,ones(1,nSamples))];



%----- espectro JRASTA: -----------

if nargout==2,

  %P=lpc2spec(a,MSE);	%--- sem hlift
  [A1,QME]=ceps2lpc(c);
  P=lpc2spec(a,QME);	%--- com liftering

end












%-- test---
%for k=1:Q,
% plot(1:72,c(k,:),1:72,c2(k,:)),pause
%end
%----------- tests: comparison with Morgan's results ----------
%rasta -J -C -w 32 -i c:\matlab\exps\dig.dig -o c:\matlab\exps\dig.cep
%c2=readFloat32('dig.cep');	%--- output
%c2=reshape(c2,Q,nSamples);

