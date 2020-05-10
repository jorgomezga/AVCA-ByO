  function [f,q,w,k]=MeddisHC(s,fs,A,B,h0,f_ini)
% function [f,q,w,k]=MeddisHC(s,fs,A,B,h0,f_ini);
%
% Meddis IHC/Synapse model
%
% s:	input signal
% fs:	sampling frequency
% A,B,h0: model constants that set fspo, fsat and fiber threshold.
% (see MedConst.m)
%
% If f_ini (mean rate output at t=0) is given,
% then use initial conditions that set f(t=0)=f0,
% else use initial conditions for sustained zero input (s=0 and f_ini=fspo)
%
% Model diagram:
%
% s(t)     z(t) +----------+ k(t)      +----+   c(t)
%   --->(+)---->|g0*z/(z+B)|---->(x)---| H1 |---+-----(x)---> f(t)
%        |      +----------+      |    +----+   |      |
%        |                   q(t) |             |      |
%                                 | -  +----+   |
%        A               M=1 --->(+)<--| H2 |<--+      h0
%                               +    r +----+
% Default constants:
% A=5; B=300; h0=50000; g0=2000;
% y0 = 5.05; l0 = 2500; r0 = 6580; x0 = 66.31; fs=20 kHz
%
% Ref:
% Meddis, "Simulation of Auditory-Neural Transduction: Further Studies", 
% JASA 83(3), pp. 1056-1063,1988. Also  JASA 87(4), pp. 1813-1816, 1990

% F. Perdigao, Coimbra 1996


if nargin<2, fs=20000; end
if nargin<3, 
  A=5;B=300; h0=50000;
  disp('Default Model!');
end

N=length(s);
T=1/fs;

%--- model constants -------

C = A+B;
g0=2000;	%--- This value is also fixed.
		%--- Consts. that may vary for desired fsat,fspo: A,B,h0
y0 = 5.05;
l0 = 2500;		% H1(z)=C(z)/E(z); e(t)=k(t)*q(t)=k(t)*(1-r(t))
r0 = 6580;		% H2(z)=R(z)/C(z)
x0 = 66.31;

H20 = l0/y0;		% H2(W=0)=H2(z=1)
H10 = 1/(l0+r0);	% H1(z=1)
y0T = y0*T;		% H1(z)=T*z^(-1) / ( 1 - a1*z^(-1) )
r0T = r0*T;		%
lr1T= 1-(l0+r0)*T;	%	1-(a1+a2)*z^(-1) +(a1*a2-x0*r0*T^2)*z^(-2)
x0T = x0*T;		% H2(z)= -------------------------------------------
x1T = 1-x0T;		%	   [1-a3*z^(-1)]*[1-a2*z^(-1)]
			%
			% a1=lr1T; a2=x1T; a3=1-y0*T
	
%fspo = h0/(C/A)/(g0*H10) + H20)	%-- fspo depends on A, B and h0.
%fsat = h0/(   2/(g0*H10) + H20)	%-- fsat depends on h0.
%----------------------------

if nargin>5, 	%--- initial conditions are given

  c = f_ini/h0;
  w = c*r0/x0;
  q = 1 - H20*c;

else		%--- use zero input initial conditions

  k0 = g0*A/C;
  c = k0*H10/(1+k0*H10*H20);	%-- or: c = k0*y0/(y0*(l0+r0)+k0*l0);
  w = c*r0/x0;
  q = 1 - H20*c;

end

%-- Non linearity -----------

kT=max(s+A,0);
kT = g0*T*kT./(kT+B);	%--- kT(n)=T*k(n)

%=============================================
if nargout ==1,	%---- The only output is f(n)

f=c;

for n=1:N-1,

  e = kT(n)*q;	
  q = q + y0T*(1-q) - e + x0T*w;
  if (q>1), q=1; disp(['q(',int2str(n),') >1']), end	%--- test! (q(t) is allways < 1)
  w = w*x1T + r0T*c;
  c = c*lr1T + e;		%-- don't change order
  f(n+1)=c;
  
end

f=h0*f;
%=======================================================
else	%------- Other outputs required!


for n=1:N-1,

  e = kT(n)*q(n);	
  q(n+1) = q(n) + y0T*(1-q(n)) - e + x0T*w(n);
  w(n+1) = w(n)*x1T + r0T*c(n);
  c(n+1) = c(n)*lr1T + e;

end

f=h0*c;	   %--- Note that: e(n)=q(n)*k(n); r(n)=1-q(n); c(n)=f(n)/h0;
k=kT/T;

end %--- nargout choice




