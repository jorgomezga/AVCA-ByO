%  Geisler's 1993 model
%  Example
%
% 1. Computs Geisler's 1993 Model Impedance
% 2. Compares impedance with a 4th order rational function impedance
% 3. Computs responses with a Long-Wave model

%-------------------------------------------------------------


clc
disp('1. Computs Geisler´s 1993 Model Impedance');
disp('=========================================');


% Parameters

gama=1.29	%-------- optimal gama

Nx=500;		%-------- Number of sections
Nf=220;		%-------- Number of freq. points

L=2.4; 		%--- Cochlear Length [cm]

x=linspace(0,L,Nx);
f=logspace(log10(100),log10(30000),Nf)';	%--- from 100 to 30kHz
s=j*2*pi*f;


%--- usual calculation:
% ZG = Geis93Z(f,x,gama);

%--- or, complete:
%[ZG,ZGmod,WQ,M]=Geis93Z(f,x,gama);

[ZG,ZGmod]=Geis93Z(f,x);



disp('2. Comparison of impedance with a 4th order rational function simplified impedance')
disp('==================================================================================')
disp(' ')
disp('Impedance viewed every 50 sections as a function of frequency.')
disp('Solid lines are for Geisler´s Impedance; dashed lines for 4th order impedance.')


semilogx(f,db(ZG(:,1:50:Nx)),f,db(ZGmod(:,1:50:Nx)),':')
xlabel('f [Hz]'), ylabel('dB'), title('ZG(f)'), pause


disp(' ');
disp('Impedance viewed every 20 freq. points as a function of cochlear place.')
disp('Solid lines are for Geisler´s Impedance; dashed lines for 4th order impedance.');

plot(x,db(ZG(1:20:Nf,:)),x,db(ZGmod(1:20:Nf,:)),':')
xlabel('x [cm]'), ylabel('dB'), title('ZG(x)'), pause

disp('');
disp('Now a more detailed view at x=3L/4...')

semilogx(f,db(ZG(:,375)),f,db(ZGmod(:,375)),':')
xlabel('f [Hz]'), ylabel('dB'), title('ZG(f), x=1,8cm'), pause



clc
disp('3. Computs responses with a Long-Wave model')
disp('===========================================')

%--------------- Transmission Line with a BOX MODEL -------------

disp(' ');
disp(' Computation of partition displacement with a BOX model')
disp(' Stapes velocity normalized...')


%uf=ones(size(f));
ux=ones(size(x));


Dx=L/Nx;		%---- sections width
rho=1;			%---- fluid density
h=0.1;			%---- efective height of scala
zs=2*s*rho*Dx;		%---- series impedance (fluid mass)
zp=h*ZG/Dx;		%---- parallel impedance
%----------------

[P,Zin]=TLineVs(1,zs,zp);	%---- stapes normalized velocity

VBM=-P./ZG;		%--- partition velocity
DBM=VBM./(s*ux);	%--- partition displacement

disp('Displacemet as a function of x for several freq.s')

%plot(x,db(P(1:40:Nf,:))),pause
%plot(x,db(VBM(1:40:Nf,:))),pause
plot(x,db(DBM(1:40:Nf,:)))
title('DBM(x)'), xlabel('x [cm]'), ylabel('dB'), pause

disp('Displacemet as a function of frequency every 50 sections')

semilogx(f,db(DBM(:,1:50:Nx)))
title('DBM(f)'), xlabel('f [Hz]'), ylabel('dB'), pause
semilogx(f,fase(DBM(:,1:50:Nx)))
title('arg{DBM}'),xlabel(' f [Hz]'),ylabel('x pi'), pause


clc
disp(' Cochlear input impedance')

semilogx(f,db(Zin))
title('Cochlear input impedance'), xlabel('f [Hz]'), ylabel('dB'), pause





%---------------POWER ABSORBED------------------------

clc
disp(' Ratio of power absorbed by partition and input power for 2, 4 and 8 kHz')

Vs=1;
Ps=Vs*zs+P(:,1);	%--- OR: Ps=Vs*Zin
%%V=([Ps,P(:,1:Nx-1)]-P)./(zs*ux);	%--- Fluid Power
Potin=real(Ps.*conj(Vs))/2;
%%semilogx(f,Potin)
PotBM =real(P.*conj(P./zp))/2;	%--- Power absorbed by partition

for k=[116,143,169];
 plot(x/L,(PotBM(k,:)/Potin(k))),hold on
end
hold off, xlabel('x/L'),ylabel('dB'),title('Pot/Potin')

%=======================================================


return


%--------- Parameter variation with gama: ---------------

Ng=100;
g=linspace(1.4,0.1,Ng);

[WZag,WHzg,WHpg,QZag,QHzg,QHpg]=Geis93G(g);

plot(g,WZag), title('WZa(gama)'),xlabel('gama'),pause
plot(g,WHzg,g,WHpg), title('WHz(gama), WHp(gama)'),xlabel('gama'), pause
plot(g,QZag,g,QHzg,g,QHpg), title('Q(gama)'),xlabel('gama'), pause



