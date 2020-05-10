  function [H,E]=PLPfilt(N,Nband,fs,flag_pe)
% function [H,E]=PLPfilt(N,Nband,fs,flag_pe)
%
% Defines the Nband filters for use in PLP analysis
% (H. Hermanski, JASA 87(4), 1990)
%
% N: number of freq. points in the range [0..pi] (pi included!! Use N=Nfft/2+1) 
% Nband: Num. of effective critical band filters (excluding CFs of 0 and pi)
% fs: sampling frequency
% H: matrix of filter responses (Nband x N)
% E: Preemphasis factors (Nband x 1). 
% If "flag_pe" is input and is different from 0, 
%    apply preemphasis factors to the filters.
%  NOTE: for RASTA, preemphasis must be computed later,
%        so E is also returned.
%      
% Defaults: Nband=16; fs=10000

% Fernando Perdigao, Dec. 97
% fp@it.uc.pt


if nargin<1, error('no. of freq. points must be specified'); end
if nargin<4, flag_pe=0; end
if nargin<3, fs=10000; end
if nargin<2, Nband=16; end


W=(0:N-1)*pi/(N-1);	%--- NOTE: W on range [0..pi]
f=fs*W/(2*pi);
w=2*pi*f;
%--- bark scale ------
bark = 6*log( w/(1200*pi) + sqrt((w/(1200*pi)).^2+1) );

%---- Nyquist freq in barks --------------
bark_fN = 6*log( fs/(1200) + sqrt((fs/(1200)).^2+1) ); % for fN=fs/2

bark_step = bark_fN/(Nband+1);
zc=(1:Nband)'*bark_step; 	%-- CFs in Bark scale

%---- CFs in linear freq. -------
fsq=(600*sinh(zc/6)).^2;
%w=1200*pi*sinh(zc/6);
%---- Equal-Loudness preemphasis ------
%E=(w.^2+56.8e6).*(w.^4)./((w.^2+6.3e6).^2)./(w.^2+0.38e9);
E=((fsq./(fsq+1.6e5)).^2).*(fsq+1.44e6)./(fsq+9.61e6);
%E=E/max(E); %-- normalized

H=zeros(Nband,N);

for k=1:Nband,
 z0=zc(k);
 z1=z0-0.5;
 z2=z0+0.5;
 zinf=z0-2.5;
 zsup=z0+1.3;
 i=find( bark>zinf & bark< z1 ); 
   H(k,i)= 10.0.^(bark(i)-z1);		%--- low freq. side
 i=find( bark>z2 & bark< zsup ); 
   H(k,i)= 10.0.^(-2.5*(bark(i)-z2));	%--- high freq. side
 i=find( bark>=z1 & bark<= z2 );
   H(k,i)= ones(size(i));
 if flag_pe, H(k,:)=H(k,:)*E(k); end
end







