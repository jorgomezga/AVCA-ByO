  function PSD=lpc2spec(a,EMQ,NPSD,fs,range)
% function PSD=lpc2spec(a,EMQ,NPSD,fs,range);
%
% Computs the LPC spectrum given a matrix, "a" with LP coefficients
% (in columns; a(1,:)=1)
% returns matrix PSD (NPSD x Nf) with LP power spectrum, where:
%   NPSD - no. of points in frequency, [0 .. pi] (default: 128)
%
% fs - sampling freq. (Opcional; if specified shows sonogram)
% range - sonogram range relative to max magnitude (default=100dB)
%
% PSD=|EMQ/A(W)|^2, where: 
%     A(z)=sum(ak*z^(-k)), k=0..p
%

[p,Nf]=size(a); p=p-1;

%------- exp(-jWn) ---------------

if nargin<3, NPSD=128; end
if nargin<5, range=100; end

W=(0:NPSD-1)'*pi/(NPSD-1);	%--- NPSD+1 freq. values in range [0..pi].
expW = exp(-j*W*(0:p));	%--- exp(-j*k*W), size(expW)=(NPSD,p+1)

%------- LPC ---------------------

PSD=zeros(NPSD,Nf);

for n=1:Nf,

 H = EMQ(n)./(expW*a(:,n));	%--- H(W)=1/sum(ak*exp(-j*k*W))= 1/A(W);
 PSD(:,n)=H.*conj(H);	%--- PSD = |H(W)|^2  (columns per frame)

end

if nargin>3,				%--- plot sonogram LPC
  f=linspace(0,fs/2,NPSD);
  t=1:Nf;
  min_P=max(max(PSD))/10^(range/20);
  imagesc(t,f,db(PSD,min_P))
  axis xy;
%  colormap(flipud(gray))
  xlabel('Time')
  ylabel('Frequency')
  grid off
end
if nargout==0,
  fprintf('No. of frames: %d\n',Nf);
  clear PSD
end

