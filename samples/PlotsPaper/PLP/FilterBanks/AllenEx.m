% AllenEx
% Computes cochlear pressure, P, cilia velocity, Vc, with Allen's model.
% see also AllenZ.m
% Ref: S. Puria, J. Allen, "A Parametric Study of Cochlear Input Impedance", 
%      JASA 89(1), pp. 287-309, Jan. 1991.


clear
clc
disp('Computes cochlear pressure, P, cilia displacement, Vc, with Allen´s model')

M=0.04;
Mb=M/1.2;
Mt=0.2*Mb;
L=2.5;


disp('Uses 1200 sections and 100 freq. points log spaced from 100Hz to 20kHz')
Nx=1200; 
Nf=100; 
L=2.5;
x=linspace(0,L,Nx);
f=logspace(log10(100),log10(20000),Nf)';


disp(' Impedance calculation...')


[Zp,Ht]=AllenZ(f,x);

s=j*2*pi*f;

disp(' Transmission Line box model with normalized stapes velocity...')

Dx=L/Nx;
h=0.1;		%--- effective height of the scalae [cm] 
rho = 1;	%--- fluid density [g/cm^3]
zs = 2*s*rho*Dx;
zp = h*Zp/Dx;

[P,Zin]=TlineVs(1,zs,zp);	%---- stapes velocity normalized 
VBM = -P./Zp;
Vc = Ht.*VBM;

disp(' Cilia velocity display for 3 positions at 1.67, 1.87 and 2.1 cm ')
semilogx(f,db(Vc(:,[800,900,1000])))
title('|Vc|'),xlabel('f [Hz]'), ylabel('dB'), pause

semilogx(f,fase(Vc(:,[800,900,1000])))
title('arg{Vc}'),xlabel('f [Hz]'), ylabel('x pi')

 
