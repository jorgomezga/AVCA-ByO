  function [Zp,Hi] = Neely93Z(f, x, gamma)
% function [Zp,Hi] = Neely93Z(f, x, gamma)
% Neely's(1993) Coclear Model. Computation of partition impedance Zp(f,x)
% Ref: JASA 94 (1) July 93 pp 137
%
% Input:
% f - column vector for required frequencies
% x - row vector for x positions (must be between 0 and L=2.5cm)
% gamma - feedback gain (default: gamma=1)
% Output:
% Zp - Model Partition SPECIFIC Impedance ((Nf x Nx) elements)
% Hi - IHC gain function

% F. Perdigao, Coimbra, Feb. 1996 (fp@it.uc.pt)
%==============================================


if nargin==2, gamma=1; end	%-- sets the default gamma

[n,Nx]=size(x);
if (n>1),
  error('--- x must be a row ---')
end
[Nf,n]=size(f);
if (n>1),
  error('--- f must be a column ---')
end

if length(gamma)>1,
 error('--- only one gamma value allowed ---')
end

%-------------------- INITIALIZATIONS -----------------------------

rho =1;		%--- cgs units
L=2.5;		%--- cochlea length  (cm)
wHB=2*pi*1390;  %--- cutoff freq. of Hair Bundle
%bw=0.01;	%--- effective width of BM
%Ap=bw*L/Nx;  	%--- partition area (=bw*delta)
gr=0.1;

%------------------------------------------------------------------

%---- coeficient variation with x (mechanical parameters) ----
%---- Parameters taken from Table I of JASA paper ------------


xx=[0,L/2,L];

Kb = exp(polyval(polyfit(xx,log([1.14e8,  4.19e6,  5.97e4 ]),2)  ,x));
Kt = exp(polyval(polyfit(xx,log([1.99e4,  2.21e4,  3.16e4 ]),2)  ,x));
K0 = exp(polyval(polyfit(xx,log([1.05e4,  9.23e3,  1.25e4 ]),2)  ,x));
Mb = exp(polyval(polyfit(xx,log([9.14e-6, 9.6e-6,  1.06e-5]),2)  ,x));
Mt = exp(polyval(polyfit(xx,log([5.64e-4, 1.02e-3, 1.06e-2]),2)  ,x));
Rb = exp(polyval(polyfit(xx,log([2.08e-2, 2.03e-2, 1.88e-2]),2)  ,x));
Rt = exp(polyval(polyfit(xx,log([149,     63.4,    27     ]),2)  ,x));
R0 = exp(polyval(polyfit(xx,log([2037,    282,     38     ]),2)  ,x));
gf = exp(polyval(polyfit(xx,log([1.42e5,  1.05e4,  3.68e2 ]),2)  ,x));
tauf=exp(polyval(polyfit(xx,log([1.4e-4,  6.92e-4, 5.29e-3]),2)  ,x));
taur=exp(polyval(polyfit(xx,log([1.35e-4, 3.61e-4, 2.5e-3 ]),2)  ,x));
%Ac= exp(polyval(polyfit(xx,log([5.52e-3, 3.17e-3, 4.27e-3]),2)  ,x));
alfa =       polyval(polyfit(xx,[0.081,   0.035,   0      ] ,2)  ,x);



w=2*pi*f;
s=j*w;
invs=1 ./s;
uf=ones(size(f));
ux=ones(size(x));

%------- Specific impedances: -------------------------

Zb=invs*Kb + uf*Rb + s*Mb;		%--Zb=(Kb./s + Rb + s*Mb );
Zt=invs*Kt + uf*Rt + s*Mt;		%--Zt=(Kt./s + Rt + s*Mt );
Z0=invs*K0 + uf*R0;			%--Z0=K0./s+R0;

Hc=gamma*gr*(uf*gf)./((1+s*tauf).*(1+s*taur));
H0 = 1.0./(1+Hc+Z0./Zt);

Zp = Zb + H0.*Z0;


%----------------------------------------------------

%-- NOTE: the partition acoustic impedance is now Zp/Ap.
%-------- the mechanical impedance is Zp*Ap

if nargout < 2, return, end;

Hi = H0.*( 1 + (uf*alfa).*Hc ).*( (s./(s+wHB))*ux );


