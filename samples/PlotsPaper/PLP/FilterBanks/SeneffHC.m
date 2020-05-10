  function [f,z,c,v,q]=SeneffHC(x,fs, A, B, GHW, f_ini)
% function [f,z,c,v,q]=SeneffHC(x,fs, A, B, GHW, f_ini);

%  Implement Seneff's hair cell model (stage II)
%
%       +----+ s +         +-----+    z(t) +-----+   v(t)      |\ AGC
% x(t)--| NL |--->(+)-(x)--| HWR |---*-----| LPF |---*---------| >----->f(t)
%       +----+     |   |   +-----+   |     +-----+   |         |/|
%                - |  wa             |               |           | 
%              c(t)|    +------+     |               |  +---+    |
%                  +----| H(z) |-----+               +--|LPF|----+
%                       +------+                        +---+  q(t)   
%
% Original results with fs=16KHz
% Forward differences are used to aprox. the derivative.
%
% A, B, GHW: parameters that define spontaneous and saturated rate
% f_ini: initial value of f(t) (for filter initial conditions)
%
% see also: SenConst.m; Slaney's Auditory toolbox.

%--- Model constants------------

if nargin==1, fs=16000; end
if nargin==2,
 A=10;
 B=65;
 GHW=2.35;
end



%--- IHC model time constants

T=1/fs;
tau1=15/1000;	%--- time constant of 15 ms
tau2=120/1000;	%--- time constant of 120 ms
wb=1/tau2;
wa=1/tau1 - wb;	%--- tau1 = 1/(wa+wb)
bH=1-wb*T;	%--- pole of H(z)=C(z)/Y(z) = [T*z^(-1)]/[1-bH*z^(-1)]
		%--- (using forward differences). Note: H(z=1)=1/wb.


%--- FILTERS: LPF and AGC -----------

lpAlpha = exp(-1/(fs*0.04e-3));		%-- tauLP=0.04ms;
aLPF = poly([lpAlpha lpAlpha lpAlpha lpAlpha]);
bLPF = sum(aLPF);

%---- AGC --------

alpha_agc = exp(-1/(fs*0.003));		%--- tau_agc = 3ms
KAGC = 0.002;
bAGC=1-alpha_agc;
aAGC=[1,-alpha_agc];


%------- Initial Conditions ------------------
if nargin==6,
 uz = f_ini/(1-f_ini*KAGC*1.01);
else
 uz = GHW*wa*wb/(wa+wb);	%--- rest conditions for null input
end

zini=filtic(bLPF,aLPF,uz*ones(size(aLPF)),uz*ones(size(bLPF)));
qini=filtic(bAGC,aAGC,uz*ones(size(aAGC)),uz*ones(size(bAGC)));
c = uz/wb;

%=========================================================
%---- Input of non-linearity -------------------------
%=========================================================


s=GHW*(A*atan(B*max(0,x))+exp(A*B*min(0,x)));

%=========================================================
%--- short-term adaptation -------------------------------
%=========================================================
N=length(x);

if nargout <2,		%--- if only y(n) and z(n) is needed...
for n=1:N
  z(n) = max(0,wa*(s(n)-c));
  c = bH*c + z(n)*T;
end
else			%--- other than y(n) ...
for n=1:N
  z(n) = max(0,wa*(s(n)-c(n)));
  c(n+1) = bH*c(n) + z(n)*T;
end
c(N+1)=[];	%--- exclude last sample

end %--- nargout choice
%=========================================================
%--- Low Pass filter (Loss of Synchrony)
%=========================================================

v = filter(bLPF,aLPF,z,zini);

%=========================================================
% Adaptive Gain Control: f(n) = v(n)/(1 + KAGC*q)
%=========================================================

q = filter(bAGC,aAGC,v,qini);
f = v./(1+KAGC*q);

