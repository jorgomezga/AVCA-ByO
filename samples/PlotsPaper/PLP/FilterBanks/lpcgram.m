  function [PSD,Nf]=lpcgram(x,p,N,M,NPSD,fs,range)
% function [PSD,Nf]=lpcgram(x,p,N,M,NPSD,fs,range)
%
% Computs and displays (order p) LPC sonogram given signal x
% Uses:
%  preemphasis (factor 0.95)
%  Hamming window of length N every M samples
%
% Returns PSD, the power spectrum, with NPSD (default=N) freq. samples in [0..pi]
% Returns also the number of frames (the number of column in matrix PSD).
%
% fs: sampling freq. If specified, displays the LPC sonogram in dB
% range: dB sonogram range from max value (default: 100dB)
%
% PSD=|EMQ/A(W)|^2,
%     A(z)=sum(ak*z^(-k)), k=0..p
%
%	Example:
%	fs=10000; x=sin(2*pi*1000/fs*(0:999));
%	PSD=lpcgram(x,10,256,64,100,fs);
%       pause, plot(db(PSD))
%
%--------------------------------------------------------------------------


if nargin<4, error('too few parameters. Minimal call: lpcgram(x,p,N,M)'); end
if min(size(x))>1, error('too few samples'); end
if nargin==4, NPSD=N; end
if nargin<7, range=100; end

nx=length(x);		%--- num. de amostras de x.
x=x(:);
x=x-0.95*[0;x(1:nx-1)];	%--- preemphasis


Nf = fix((nx-N)/M)+1	%--- no. of frames

%--- block signal into frames

y =zeros(N,Nf);
ind=ones(N,Nf);
ind(1,:)=(0:Nf-1)*M+1;	%-- frame index begining: 1;M+1;2*M+1;...
ind=cumsum(ind);	%-- indexes from x.
y(:)=x(ind);		%-- y: matrix where each column is a signal frame.

h = 0.54 - 0.46*cos(2*pi*(0:N-1)'/(N-1));	%-- Hamming window
y=y.*h(:,ones(1,Nf));	%-- windowing. 


%------ for computation of autocorrelation vector --------

iR=cumsum([(1:N);ones(p,N)]);	%--- indexes for autocorrelation R
R=zeros(size(iR));		%--- shape R
zp=zeros(p,1);


%------- exp(-jWn) ---------------

W=(0:NPSD-1)'*pi/(NPSD-1);	%--- NPSD freq. values in [0..pi].
expW = exp(-j*W*(0:p));		%--- exp(-j*k*W), size(expW)=(NPSD,p+1)

%------- LPC ---------------------

PSD=zeros(NPSD,Nf);

for n=1:Nf,

 z=[y(:,n);zp];
 R(:)=z(iR);
 r=R*y(:,n);		%---- autocorrelation

 a = [ 1; -toeplitz(r(1:p))\r(2:p+1)];	%--- LPC vector (column)
 EMQ=r'*a;				%--- quadratic mean error

 H = sqrt(EMQ)./(expW*a);	%--- H(W)=1/sum(ak*exp(-j*k*W))= 1/A(W);
 PSD(:,n)=H.*conj(H);		%--- PSD = |H(W)|^2  (colunas por frame)

end

%--- plot LPC sonogram

if nargin>5,
  f=linspace(0,fs/2,NPSD);
  t=linspace(0,nx/fs,Nf);
  min_P=max(max(PSD))/10^(range/20);
  imagesc(t,f,db(PSD,min_P))
  axis xy;
%  colormap(flipud(gray))	%--- gray scale
  colormap(jet)
  xlabel('Time')
  ylabel('Frequency')
  grid off
end
if nargout==0,
  clear PSD Nf
end

