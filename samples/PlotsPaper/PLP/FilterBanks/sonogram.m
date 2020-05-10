  function sonogram(x,N,M,Nfft,fs,range)
% function sonogram(x,N,M,Nfft,fs,range)
%
% Sonogram of signal x, computed with N sample Hamming window every M samples 
% and FFT power spectrum, with Nfft/2 freq. points
% Inputs:
%   x : signal (vector with Nx samples)
%   N : window length
%   M : window lead (frame rate in samples)
%
%   Nfft: No. of FFT points, [0..2*pi[ (default: pow2(nextpow2(N)))
%   fs: sampling frequency (optional)
%   range : sonogram range (in dB) for ploting (default 80 dB)
%
% see also lpcgram.m


if nargin<3, error('x,N or M not specified')
elseif nargin<4,
   Nfft=pow2(nextpow2(N));		%--- default no. of FFT points
   range=80;
elseif nargin<5,
   Nfft=pow2(nextpow2(Nfft));		%--- no. of FFT points
   range=80;
elseif nargin<6,
   range=80;
end;
if min(size(x))~=1, error('x must be a vector'); end



nx=length(x);
Nf = fix((nx-N)/M)+1;	%--- no. of frames
disp(['No. of frames: ',int2str(Nf)]);


y =zeros(N,Nf);
ind=ones(N,Nf);
ind(1,:)=(0:Nf-1)*M+1;	%--- frame starting index
ind=cumsum(ind);	%-- ind: x indexes for frames

y(:)=x(ind);		%-- frame matrix (each column is a frame)
h = 0.54 - 0.46*cos(2*pi*(0:N-1)'/(N-1));	%-- Hamming window
y=y.*h(:,ones(1,Nf));						%--- Apply hamming window

y(N+1:Nfft,:)=zeros(Nfft-N,Nf);	%--- zero-pad

y=fft(y);
y=y(1:(Nfft/2),:);		%--- for W=[0..pi[
y=y.*conj(y);		%--- power spectrum, |X(W)|^2


  range=10^(range/20);
  min_y = max(max(y))/range;
  if ~isempty(fs),
   f=(0:Nfft/2-1)*fs/Nfft;	%--- [0..fs/2[
   t=((N:M:nx)-N/2)/fs;		%--- time at window center
   imagesc(t,f,db(y,min_y))
   xlabel('Time [s]')
   ylabel('Frequency [Hz]')
  else
    imagesc(db(y,min_y));
  end
 grid off
 axis xy;

% colormap(flipud(gray))
  colormap(jet)
