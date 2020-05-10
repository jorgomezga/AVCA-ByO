%---- Kanis-de Boer's model (linear active case)
%
% Computes Pressures and Velocities for the Kanis & de Boer' model
%
% See also: KansiZ.m  KanisG.m
%-------------------------------------------------------------------------------------

clear
clc
disp(' Kanis and de Boer´s model')
disp(' Takes 500 sections and 220 freq. points from 200Hz to 20kHz')

Nx=500;
Nf=220;
L=0.032;	%--- coclhear length in [m]

x=linspace(0,L,Nx);
f=logspace(log10(200),log10(20000),Nf)';
Zp=KanisZ(f,x);
%--------------------------------------------------

disp('Transmission Line Box model with normalized stapes impedance')

Dx=L/Nx;
h=0.001;	%--- MKS units
rho = 1000;
s=j*2*pi*f;

zs = s*rho*Dx;
zp = h*Zp/(2*Dx);

[P,Zin]=TlineVs(1,zs,zp);

VBM=-2*P./Zp;

%-----------------------------------------------------------

disp('Basilar Membrane (BM) Velocity as a function of x for several freq.s')

plot(x,db(VBM(1:40:Nf,:)))
title('VBM(x)'), xlabel('x [m]'), ylabel('dB'), pause

disp('BM velocity as a function of frequency every 50 sections')

semilogx(f,db(VBM(:,1:50:Nx)))
title('VBM(f)'), xlabel('f [Hz]'), ylabel('dB'), pause
semilogx(f,fase(VBM(:,1:50:Nx)))
title('arg{VBM}'),xlabel(' f [Hz]'),ylabel('x pi'), pause


clc
disp(' Cochlear input impedance')

semilogx(f,db(Zin))
title('Cochlear input impedance'), xlabel('f [Hz]'), ylabel('dB'), pause



%---------------POWER------------------------
clc
disp(' Ratio of power absorbed by partition and input power for 2, 4 and 8 kHz')


Vs=1;
Ps=Vs*Zin;
Potin=real(Ps.*conj(Vs))/2;
semilogx(f,Potin)
PotBM =real(P.*conj(P./zp))/2;	%--- Power absorbed by partition

for k=[110,143,176],
 plot(x/L,(PotBM(k,:)/Potin(k))),hold on
end
hold off, xlabel('x/L'),ylabel('dB'),title('Pot/Potin')

return

%==================== variation with gama ==============================


[ZKB,WQ,M]=KanisZ(f,x);		%-- ZKB: partition impedances; M: mass (constant)

%--- test ----
semilogx(f,db(Zp(:,100:100:Nx)),f,db(ZKB(:,100:100:Nx)),':')


%------- WQ = [WZa,WHz,WHp,QZa,QHz,QHp]; --------
WZa=WQ(1); WHz=WQ(2); WHp=WQ(3);
QZa=WQ(4); QHz=WQ(5); QHp=WQ(6);


gama=linspace(0.01,1,50)';
[WZag,WHzg,WHpg,QZag,QHzg,QHpg]=KanisG(gama,L/2);


plot(x/L,log10(WZa*w0),x/L,log10(WHz*w0),x/L,log10(WHp*w0),gama,log10(WZag*w0(Nx/2)))



%---- Assimptotic Zin -------------

S0 = 1e10;	%-- stiffness value at x=0
alfa = 300;	%-- m^-1; MKS units

Rin = sqrt(rho*h*S0/2);
v0=(2/alfa)*(2*pi*f)*sqrt(2*rho/h/S0);
Zass=j*Rin*(BesselJ(0,v0)-j*BesselY(0,v0))./(BesselJ(1,v0)-j*BesselY(1,v0));

semilogx(f,db(Zin),f,db(Zass))



