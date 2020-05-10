  function [f0,x]=NeelyCF(Nx,gamma,Zh)
% function [f0,x]=NeelyCF(Nx,gamma,Zh)
%
% Computs the characteristic frequencies in the Neely's 1993 model
% as a function of x, the cochlear position
% The characteristic frequencies are defined as the frequencies for which
% the IHC displacement responses peak (Best Frequencies).
%
% Input vars:
%	Nx (Number of cochlear sections. Default: 500)
%	gamma (default 1)
%	Zh (optional, the helicotrema impedance (for no reflexions))
%
% NOTE1: Uses 220 frequency values from 100Hz to 30000 Hz.
% NOTE2: Uses polynomial interpolation to smoth the f0 curve
%
% To see CF´s as a function of x, use: [f0,x]=NeelyCF; plot(x,f0)
% To compare f0 with Liberman function (cat data), use:
% fcat=456*(10 .^(2.1*(1-x/L)) -0.8); plot(x,log10(f0),x,log10(fcat),':')



L=2.5;
Nf=220;
fmin=100;
fmax=30000;
f=logspace(log10(fmin),log10(fmax),Nf)';
s=j*2*pi*f;


if nargin==0,
 Nx=500;
 gamma=1;
elseif nargin<2, gamma=1; end;
end

x=linspace(0,L,Nx);


%--- Partition Impedances and IHC transfer function Hi: -----------------

[Zp,Hi] = Neely93Z(f,x,gamma);
size(Zp)
size(Hi)
pause


%-------------- middle-ear impedance -------------------------------
As=1e-2; 
Ae=15e-2; 
gm=0.5; 
Mm=5e-4; % <----- This sould be the correct value for Mm
Rm=15; 
Km=1.5e5;

Zms = (Km./s + Rm + s*Mm)/((gm*As)^2);	%--- Zm reduced to transformer secundary
Pes = Ae/(gm*As);			%--- with  Pe=1 (pressure at eardrum)
%--------------------------------------------------------------------------

%--- impedance at helicotrema for small reflections -----------------------

if nargin<3,
 Zh = N93ZL(f,gamma);	
end

%---------- series and parallel impedances for Long Wave model ------------
Dx=L/Nx;
bw=0.01;	%--- effective width of BM
Ap=bw*Dx;  	%--- partition area
Rf=2e-5;	%--- There must be a sign error in the paper
rho=1;
Ac= exp(polyval(polyfit([0,L/2,L],log([5.52e-3, 3.17e-3, 4.27e-3]),2)  ,x));
Mf=2*rho*Dx./Ac;

size(Zp)
zs = Rf + s*Mf;			%--- Rf does little or nothing
zp = [Zp(:,1:Nx-1)/Ap,Zh];	%--- Zh: for small reflections at helicotrema


P=TlineME(Pes,Zms,zs,zp);	%--- Transmission Line Model


%--------- BM velocity and IHC displacement -------------------------------------

Vb=-P./Zp;			%--- eq. (5) and (12) in Neely's paper
Di = Hi.*Vb./s(:,ones(1,Nx));	%-- displacement of IHC (eq. 24)

%--------------------------------------------------------------------------------

%--- test --------
%semilogx(f,db(Di(:,1:50:Nx))), title('Displacement of IHC'), xlabel('f[Hz]'),pause


%----- Coclear map: ------------------------------------------------------------

for k=1:Nx,
 aux=db(Di(:,k));
 i=find(aux==max(aux));
 f0(k)=f(i(1));
end

%---- Polynomial interpolation: ---------------

i=find(x>0.34);
c=polyfit(x(i),log10(f0(i)),4);
%plot(x,log10(f0),x,polyval(c,x)),pause
%plot(x,log10(f0)-polyval(c,x)),pause
f0=10 .^(polyval(c,x));

%----------------------------------------------------------

disp('To see the CF´s as a function of x, use:')
disp('[f0,x]=NeelyCF; plot(x,f0)')
disp(' ')
disp('To compare f0 with Liberman function (cat data), use:')
disp('fcat=456*(10 .^(2.1*(1-x/2.5)) -0.8);')
disp('plot(x,log10(f0),x,log10(fcat),'':'')')



