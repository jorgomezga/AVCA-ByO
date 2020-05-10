% Neely and Kim's (1986) cochlear model 
% Computes Pressures and Velocities
% Ref: S. Neely, D. Kim, "A Model for Active Elements in Cochlear Biomechanics", 
%J. Acoust. Soc. Am. 79(5), pp. 1472-1480, May 1986.


clear
clc

disp(' Neely and Kim´s model (1986)');
disp(' Impedance with Nf=220 freq. points between 100 Hz and 50kHz and Nx=500 sections.');

Nf=220;
Nx=500;
L=2.5;		%--- cochlear length

f=logspace(2,log10(50000),Nf)';
x=linspace(0,L,Nx);
s=j*2*pi*f;

%[Zp,Zpmod,fQ,Mef]=NKim86Z(f,x);	%--- default gama=1
%fZa=fQ(1,:);fHz=fQ(2,:);fHp=fQ(3,:);
%QZa=fQ(4,:);QHz=fQ(5,:);QHp=fQ(6,:);
%plot(x,log(fZa),x,log(fHz),x,log(fHp))

%or simply: Zp=NKim86Z(f,x);
%--------------------------

%--- discrete computations ------
%================================

gama=1;
g=1;
b=0.4;
k1=1.1e9*exp(-4*x);
c1=20+1500*exp(-2*x);
m1=3e-3;
k2=7e6*exp(-4.4*x);
c2=10*exp(-2.2*x);
m2=0.5e-3;
k3=1e7*exp(-4*x);
c3=2*exp(-0.8*x);
k4=6.15e8*exp(-4*x);
c4=1040*exp(-2*x);


invs=1 ./s;
umf=ones(Nf,1);
umx=ones(1,Nx);

Z1=invs*k1 + c1(umf,:) + m1*s(:,umx);   %--- Z1=k1/s+c1+s*m1;
Z2=invs*k2 + c2(umf,:) + m2*s(:,umx);	%--- Z2=k2/s+c2+s*m2;
Z3=invs*k3 + c3(umf,:);			%--- Z3=k3/s+c3;
Z4=invs*k4 + c4(umf,:);			%--- Z4=k4/s+c4;
Hc=Z2./(Z2+Z3);

Zp=(g/b)*( Z1 + Hc.*(Z3-gama*Z4) );

%---------ZBM=(g/b)*Z1;
%---------ZOC=(g/b)*Hc.*(Z3-gama*Z4);



%======================================================================

disp(' Middle-ear connection')

%-------- middle ear -------------
km=2.1e5;
cm=400;
mm=45e-3;
As=0.01;
Am=0.35;
gm=0.5;
Zm = s*mm + cm + km./s;
Pi = Am/gm/As;
%---------------------------------


disp(' Transmission Line Model')

Dx=L/Nx;
h=0.1;
rho=1;
zs=2*s*rho*Dx;
zp=h*Zp/Dx;
zp(:,Nx)=sqrt(2*s*rho*h.*Zp(:,Nx))-zs;

[P,Ps,Zin]=TlineME(Pi,Zm,zs,zp);
VBM = -P./Zp;
VC=Hc.*VBM;


disp(' Ploting Responses every 100 sections...')

semilogx(f,db(P(:,1:100:Nx),1e-5))
title('|P(f)| (pressure)'), xlabel('f [Hz]'), ylabel('dB'), pause
semilogx(f,fase(P(:,1:100:Nx)))
title('arg{P(f)}'), xlabel('f [Hz]'), ylabel('x pi'), pause



semilogx(f,db(VBM(:,1:100:Nx)))
title('|VBM(f)| (BM velocity)'), xlabel('f [Hz]'), ylabel('dB'), pause
semilogx(f,fase(VBM(:,1:100:Nx)))
title('arg{VBM(f)}'), xlabel('f [Hz]'), ylabel('x pi'), pause

semilogx(f,db(VC(:,1:100:Nx)))
title('|VC(f)| (cilia velocity)'), xlabel('f [Hz]'), ylabel('dB'), pause
semilogx(f,fase(VC(:,1:100:Nx)))
title('arg{VC(f)}'), xlabel('f [Hz]'), ylabel('x pi'), pause


%---------------POWER------------------------
clc
disp(' Ratio of power absorbed by partition and input power for 2, 4 and 8 kHz')


Potin=real(Ps.*conj(Ps./Zin))/2;
PotBM =real(P.*conj(P./zp))/2;

for k=[155,131,107];
 plot(x/L,PotBM(k,:)/Potin(k)),hold on
end
hold off,  xlabel('x/L'),ylabel('dB'),title('Pot/Potin')

%=======================================================

