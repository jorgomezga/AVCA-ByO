  function [y1,y2,y3]=PLP(x,N,M,Nband,fs,p,Q)
% PLP spectral analysis of vector x
%
% P=PLP(x,N,M,Nband,fs);
%
% Performs PLP analysis on vector x, ...
% with frames of N samples, every M samples, with Nband filters.
% (H. Hermanski, JASA 87(4), 1990)
%
% x: signal
% N: frame length
% M: frame rate in samples
% Nband: num. of (effective) critical bands
% fs: sampling frequency
% P: auditory power spectrum (Nband+2)x(Nframes) : 
%    (endpoints corresponding to 0 and pi are copies from adjacent ones)
%
% Other forms:
% [a,P]=PLP(x,N,M,Nband,fs,p);
%  computs vector "a" (p+1)x(nSamples): LP coeficients (order p; a(1)=1).
%
% [c,a,P]=PLP(x,N,M,Nband,fs,p,Q);
%  computs vector "c" with Q cepstral coeficients (including c0, no liftering).  
%--- NOTE: N,M,Nband,p,Q must be integers!

% Fernando Perdigao, Dec. 97
% fp@it.uc.pt


if nargin<5, error('too few parameters...'); end

[X,nSamples,Nfft]=frames(x,N,M); %--- block into frames (with hamming window)

X=abs(fft(X));
X=X(1:Nfft/2+1,:).^2;

H=PLPfilt(Nfft/2+1,Nband,fs,1); %--- PLP filters with preemphasis


P=(H*X).^(1/3);	 %--- power spectrum with cubit root (power law)

y1=[P(1,:);P;P(Nband,:)];%--- PLP power spectrum with copies from nearest freqs.

if nargin>5,	%-- LP or cepstrum

 if p>Nband, error('LP order must be less then Nband'); end

 %--- autocorrelation: ------
 %R=real(ifft([y1;flipud(y1(2:Nband+1,:))]));
 R=real(ifft([y1;flipud(P)]));

 a=zeros(p+1,nSamples);
 MSE=zeros(1,nSamples);

  for k=1:nSamples,
    r=R(1:p+1,k);
    a(:,k) = [ 1; -toeplitz(r(1:p))\r(2:p+1)];	%-- LP coeficients
    MSE(k)=r'*a(:,k);				%--- MS error
  end

  if nargin==7,
   y3=y1;y2=a;
   y1=[log(MSE);lpc2ceps(a,Q-1)];
  else y2=y1; y1=a;
  end

end
