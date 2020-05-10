  function [y,Nf,Nfft]=frames(x,N,M,h)
% function [y,Nf,Nfft]=frames(x,N,M,h);
%
% Block a signal x into Nf frames of length N every M samples. 
% Apply a window h (1 column). Default: hamming window.
% Output:
%  matrix "y" (Nfft)x(Nf) (each column is a frame), where:
%   Nfft: next power 2 of N: Nfft=2^log2(N);
%   Nf: number of frames (or columns): Nf = fix((length(x)-N)/M)+1;

% F. Perdigao - Coimbra - Aug. 96
% fp@it.uc.pt

if nargin<3, error('too few parameters...'); end

log2N=ceil(log(N)/log(2));
Nfft=2^log2N;			%--- num. pontos p/ FFT potência de 2

nx=length(x);
Nf = fix((nx-N)/M)+1;	%--- no. of frames

y =zeros(N,Nf);		%--- signal frames by columns
ind=ones(N,Nf);
ind(1,:)=(0:Nf-1)*M+1;	%--- indexs of frame begining
ind=cumsum(ind);	%--- x indexs for frames

y(:)=x(ind);		%-- (see help colon) On the left side of an assignment statement,
			%-- A(:) fills A, preserving its shape from before.

y(N+1:Nfft,:)=zeros(Nfft-N,Nf);	%--- zero-pad

if nargin<4,
  h = 0.54 - 0.46*cos(2*pi*(0:Nfft-1)'/(Nfft-1));	%-- hamming window
end

y=y.*h(:,ones(1,Nf));	%--- apply window

