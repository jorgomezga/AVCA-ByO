  function [f,fr,v,q]=MartHC(x,fs,Th,fspo,fsat,n, f_ini)
% function [f,fr,v,q]=MartHC(x,fs,Th,fspo,fsat,n, f_ini);
%
% Martens-Immerseel IHC/Synapse Model
%
% x: input signal
% fs: sampling frequency
% Th: fiber threshold in dB; A=10^(Th/20); (default: A=max(x)/1000)
% fspo - spontaneous firing rate	(default: fspo=50)
% fsat - saturation firing rate		(default: fsat=150)
% n - power in AGC			(default: n=2)
% f_ini - Initial rate (for filter inicial conditions)
% 
% f: "instantaneous" firing rate
% fr: average firing rate
%
% Operation:
%
% v(t)=max(0,x(t)+A)
%
%           fsat*v(t)
% f(t)= ------------------
%       {B + q(t)^(1/n)}^n
%                            |\
%                            | \
% x+A   +-----+  v(t)        |  \ f(t) +-----+
% ----->| HWR |--+---------->|AGC>-----| LPF |---> fr(t)
%       +-----+  |           |  /      +-----+
%                |           | /
%                |           |/|
%                |   +-----+   |
%                +-->| LPF |---+ q(t)
%                    +-----+
% Ref:
% J. Martens, L. Immerseel, 
% "An Auditory Model Based on the Analysis of the Envelope Patterns", 
% ICASSP-90, pp. 401-404, 1990. Also:  JASA 91(6), pp. 3511-3526, Jun. 1992.

% F. Perdigao (fp@it.uc.pt)
% Coimbra, 1996




%---- defaults ----------
if nargin<3, 
 Th=db(max(x))-60; fsat=150; fspo=50; n=2;
end

tau1=0.008; %-------- 8ms
tau2=0.040; %-------- 40ms
r0 = 0.86;
c1 = exp(-1/(fs*tau1));     %---- LPF const.s
c2 = exp(-1/(fs*tau2));
bLP = [r0*(1-c1)+(1-r0)*(1-c2) , -r0*(1-c1)*c2-(1-r0)*(1-c2)*c1];
aLP = [1,-(c1+c2),c1*c2];
[bEEF,aEEF] = butter(3, 2*250/fs);	%--- Butterworth filter with fc=250Hz

A=10^(Th/20);
B = (A^(1/n))*( (fsat/fspo)^(1/n) - 1 );

%---- Rest conditions for the filters --------
if nargin <7,
 f_ini=fspo;
 q_ini=A;
else
 q_ini = ( B/((fsat/f_ini)^(1/n)-1) )^n;
end

iniq=filtic(bLP,aLP,q_ini*ones(size(aLP)),q_ini*ones(size(bLP)));
iniz=filtic(bEEF,aEEF,f_ini*ones(size(aEEF)),f_ini*ones(size(bEEF)));

v=max(0,x+A);			%---- Half-Wave Rectifier
q = filter(bLP,aLP,v,iniq);	%---- LPF for AGC

if n==1,
 f = fsat*v./(B+q);
else
 f = fsat*v./( (B + q.^(1/n)).^n );
end

if nargout>1,
 fr = filter(bEEF,aEEF,f,iniz);		%--- average firing rate (loss of synchrony)
end

