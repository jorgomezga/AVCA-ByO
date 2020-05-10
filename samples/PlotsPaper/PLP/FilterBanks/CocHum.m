  function [P,D,x,w0]=CocHum(f)
% function [P,D,x,w0]=CocHum(f)
%
% Cochlear model for the human coclea.
% Uses:
% parametric partition impedance.
% Rosowski external and middle ear model (see rosowski.m)
% Area of cochlear sections as an interpolated polynomial from
% measured data (Weber'49 in Puria and Allen, JASA 89(1),pp.287-309, 1991.)
%
% P: Pressure across cochlear partition
% D: IHC cilia displacement
% x: Coclear distance (500 points from 0 to L=3.5 cm)
% w0: Cochlear map (Liberman function)
%
% Ref:
% "Modelo computacional da coclea humana", 
% F. Perdigao, L. Sa', Acustica'98, Lisboa, 1998


function [P,D,x,w0]=CocHum(f)

[n,m]=size(f);
if min(n,m)>1, error('f must be a vector'); end
if n==1, f=f'; end
Nf=length(f);

%--- cgs units --------
rho=1;		%--- fluid mass density
M=0.1;		%--- basilar membrane mass
L=3.5;		%--- cochlear length

Nx=500;		%--- number of cochlear sections
Dx=L/(Nx);
x=linspace(0,L,Nx);

A=165.4*2*pi; a=2.1; B=0.8;
w0=A*(10 .^(a*(1-x/L))-B);	%-- Liberman function (map)


uf=ones(size(f));
ux=ones(size(x));
s=j*2*pi*f;
s2=-(2*pi*f).^2;


%---- PARAMETRIC PARTITION IMPEDANCE --------
%---- more or less according Neely's model --

waz=w0.*(0.9+1.1*x/L);
wap=0.95*w0;
daz=0.4 + 0.6*x/L;
dap=0.3 + 0.2*x/L;
QZa = exp(polyval(polyfit([0,1,L],log([15,  18,  6]),2)  ,x));
dZa=1.0./QZa;

Za=M*( s*ux + uf*(dZa.*w0) + (1.0./s)*(w0.^2) );
Ha=( s2*ux + s*(daz.*waz) + uf*(waz.^2))./( s2*ux + s*(dap.*wap) + uf*(wap.^2));
Zp=Za.*Ha;

%---- Transfer function to obtain IHC cilia velocity/displacement
%=================================================================

Hi=(s2*ux + s*(0.5*0.5*w0) + uf*((0.5*w0).^2))./( s2*ux + s*(dap.*wap) + uf*(wap.^2));



%-------- Area function ------
%-------- from Wever 1949 (see Puria & Allen)

xA=[...
0,0.05;
0.05, 0.06;
0.09, 0.064;
0.23, 0.07;
0.32, 0.056;
0.35, 0.015;
0.38, 0.014;
0.47, 0.01;
0.61, 0.0098;
0.8,  0.009;
0.95, 0.013;
1.6,  0.014;
1.7,  0.015;
1.85, 0.014;
2.4,  0.0155;
2.5,  0.012;
2.57, 0.0115;
3.02, 0.0097;
3.4,  0.008;
3.5,  0.0055];
Aorig= interp1(xA(:,1),xA(:,2),x)';

%---- polynomial inperpolation 
c=polyfit([0,0.5,1,1.5,2,2.4,3,3.5],Aorig([Dx,0.5,1,1.5,2,2.4,3,3.5]/Dx),5);
A=polyval(c,x);


%---------------------------

bw=0.1;		%--- partition width
Ap=bw*Dx;	%--- Area of partition's section
Rf=2e2;		%--- simulates fluid viscous loss


zs = Rf+ (2*s*rho*Dx)*(1.0./A);		%--- fluid acustic impedance

%---------- EXTERNAL/MIDDLE EAR (Rosowski model) ---------------

b =[1.0683e+000 -7.1108e+005  4.8589e+010 -1.9088e+015 -4.6478e+019  3.2648e+024 7.8505e+027];
a =[1.0000e+000  6.5762e+004  6.7790e+009  2.8729e+014  6.2436e+018  1.0478e+023 5.0628e+026];
Peq=freqs(b,a,2*pi*f);
Zeq1=(1e-5)*(s*1.4e6+ 2e10 + 2.3e14./s);   %---- 1e-5: passagem de MKS para cgs
b2 =[2.8939e+005,4.6250e+011,3.4374e+016,3.2750e+021,1.4567e+026,4.2793e+030,...
     8.7703e+034,1.2144e+039,9.6269e+042,4.6572e+046,5.1613e+049];
a2 =[1.0000e+000,3.3878e+005,2.3356e+010,1.9967e+015,7.2365e+019,1.6808e+024,...
     2.6765e+028,1.7384e+032,2.1299e+035,0];
Zeq=1e-5*freqs(b2,a2,2*pi*f);


%=============================================

%----- TRANSMISSION LINE model ---------

P=TlineME(Peq,Zeq,zs,Zp/Ap);
D=-P./Zp.*Hi;			%--- cilia displacement

