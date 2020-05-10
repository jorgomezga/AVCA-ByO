% Neely's Model (1993)
% Computation of Responses
% Ref: S. Neely, "A Model of Cochlear Mechanics with Outer Hair Cell Motility", 
%      J. Acoust. Soc. Am. 94(1), pp. 137-146, July 1993.


clear
clc
disp(' Neely´s Model (1993)')
disp(' Impedance computation for 220 freq. points (100Hz to 50kHz) and 500 sections')




%------------ Parameters ------------------------
Nx=500;
Nf=220;
fmin=100;
fmax=50e3;
gamma=1;	%--- feedback factor (default value)
L=2.5;

x=linspace(0,L,Nx);
f=logspace(log10(fmin),log10(fmax),Nf)';	%-- f em colunas
s=j*2*pi*f;
ux=ones(size(x));
%uf=ones(size(f));


bw=0.01;
Dx=L/Nx;
Ap=bw*Dx;
rho=1;
%---------------------------------------------------------------


[Zp,Hi] = Neely93Z(f,x);


%------------------------------------------------------------------------

%------- Simplified impedance parameters ---------------------------------

%[wZa,wHz,wHp,aZa,aHz,aHp,wz1,wz2,wp1,wp2]=N93roots(x);
%Mb = Ap*exp(polyval(polyfit([0,L/2,L],log([9.14e-6, 9.6e-6,  1.06e-5]),2)  ,x));
%Mef=Mb.*wz1./wp1;
%invs=1.0./s;
%Za = (invs)*Mef.*((s.^2)*ux + s*aZa + uf*(wZa.^2) )/Ap;
%Ha = ((s.^2)*ux + s*aHz + uf*(wHz.^2))./((s.^2)*ux + s*aHp + uf*(wHp.^2));
%Zpalt = Za.*Ha;
%semilogx(f,db(Zp(:,1:100:Nx)),f,db(Zpalt(:,1:100:Nx)),':'), pause
%clear Za Ha


disp(' Middle-ear impedance')

%------ middle-ear impedance --------------
As=1e-2; 
Ae=15e-2; 
gm=0.5; 
Mm=5e-4; % <----- This sould be the correct value for Mm
Rm=15; 
Km=1.5e5;

Zm=(Km./s + Rm + s*Mm);
Zm = Zm/((gm*As)^2);	%--- Zm reduced to transformer secundary
Pin=Ae/(gm*As);		%--- with  Pe=1 (pressure at heardrum)
%-----------------------------------------------------------------

disp(' Helicotrema impedance computation for very small reflections...')
disp(' Needs a lot of memory. Please wait')

if exist('N93ZL.mat')
  load N93ZL
else 
  ZL=N93ZL(f,gamma,Ap);
  %save N93ZL ZL
end


Rf=2e2;	%--- What is really the correct value?
Ac= exp(polyval(polyfit([0,L/2,L],log([5.52e-3, 3.17e-3, 4.27e-3]),2)  ,x));
Mf=2*rho*Dx./Ac;	%--- fluid acustic impedance

%-------------------------------

zs = Rf + s*Mf;			%--- Rf does little or nothing
zp = [Zp(:,1:Nx-1)/Ap,ZL];	%--- includes ZL for small reflections at helicotrema

%---------------------------------------------------------------------------


disp(' Transmission Line Model')

[P,Ps,ZC]=TlineME(Pin,Zm,zs,zp);	%--- Transmission Line Model

	%--- Ps: pressure at stapes
	%--- Zc: input impedance of the cochlea


Vb=-P./Zp;		%--- eq. (5) and (12) in Neely's paper
Di = Hi.*Vb./(s*ux);	%-- displacement of IHC (eq. 24)


disp(' Plot of responses (every 100 sections)')


semilogx(f,db(Vb(:,1:100:Nx)))
title('|Vb(f)| (BM velocity)'), xlabel('f [Hz]'), ylabel('dB'), pause
semilogx(f,fase(Vb(:,1:100:Nx)))
title('arg{Vb(f)}'), xlabel('f [Hz]'), ylabel('x pi'), pause

semilogx(f,db(Di(:,1:100:Nx)))
title('|Di(f)| (IHC displacement)'), xlabel('f [Hz]'), ylabel('dB'), pause
semilogx(f,fase(Di(:,1:100:Nx)))
title('arg{Di(f)}'), xlabel('f [Hz]'), ylabel('x pi'), pause


disp(' Figure 8 in Neely´s paper')

semilogx(f,-db(Di(:,[60,130,199,267,329,386,434,472]))-78)
axis([100,50000,0,110]),title('fig.8'), pause


%---------------POWER------------------------

disp(' Ratio of power absorbed by partition and input power for 2, 4 and 8 kHz')


Potin=real(Ps.*conj(Ps./ZC))/2;
PotBM =real(P.*conj(P./zp))/2;

for k=[155,131,107];
 plot(x/L,PotBM(k,:)/Potin(k)),hold on
end
hold off,  xlabel('x/L'),ylabel('dB'),title('Pot/Potin')
%========================================================


