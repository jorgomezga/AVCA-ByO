  function [Di,f,x,Zc]=NeelyDi(Nf,fmin,fmax,gamma,Zh)
% function [Di,f,x,Zc]=NeelyDi(Nf,fmin,fmax,gamma,Zh);
%
% Computs IHC displacement frequency curves for the Neely's cochlear model -------
% Returns also the input impedance of the cochlea
% Input: gamma (0.2-1.0) (default 1)
%	 Nf - no. of frequency points (default 300)
%	 fmin - lowest freq. (default 100 Hz)
%	 fmax - highest freq. (default 30000 Hz)
%	 Zh - helicotrema impedance (optional)
%
% Output:  Di - IHC displacement
%          f  - vector of Nf freq. points (log distributed)
%	   x  - vector of 500 points between 0-L (linear scale)	   
%	   Zc - Cochlear input impedance
%
% Example: [Di,f]=NeelyDi(220); semilogx(f,db(Di(:,1:50:500)))
% see also Neely93Z TLineME Neely93
%

%-- Fernando Perdigao - DEEUC/IT
%-- 1st version: Aug.96


if nargin==0,
  Nf=300;
  fmin=100;
  fmax=30000;
  gamma=1;
elseif nargin==1,
  fmin=100;
  fmax=30000;
  gamma=1;
elseif nargin==2,
  fmax=30000;
  gamma=1;
elseif nargin==3,
  gamma=1;
end

%--------------------- model costants ----------------------
Nx=500;
L=2.5;
As=1e-2; 
Ae=15e-2; 
gm=0.5; 
Mm=5e-4; % <----- This sould be the correct value for Mm
Rm=15; 
Km=1.5e5;

Dx=L/Nx;
bw=0.01;	%--- effective width of BM
Ap=bw*Dx;  	%--- partition area
Rf=2e-5;	%--- There must be a sign error in the paper
rho=1;

%-------------------------------------------------------------

x=linspace(0,L,Nx);

f=logspace(log10(fmin),log10(fmax),Nf)';
s=j*2*pi*f;

if nargin <5,
 disp('Computing Zh ...') 
 Zh = N93ZL(f,gamma);
 disp('Zh calculated.')
end

[Zp,Hi] = Neely93Z(f,x,gamma);


%-------------- middle-ear impedance --------------

Zms = (Km./s + Rm + s*Mm)/((gm*As)^2);	%--- Zm reduced to transformer secundary
Pes = Ae/(gm*As);			%--- with  Pe=1 (pressure at eardrum)



%--------------- section's impedances -----------------------

Ac= exp(polyval(polyfit([0,L/2,L],log([5.52e-3, 3.17e-3, 4.27e-3]),2)  ,x));
Mf=2*rho*Dx./Ac;

zs = Rf + s*Mf;			%--- Rf does little or nothing
zp = [Zp(:,1:Nx-1)/Ap,Zh];	%--- Zh for small reflections at helicotrema

%---------------- Transmission line model -------------------

[P,Ps,Zc]=TlineME(Pes,Zms,zs,zp);	%--- Transmission Line (Long Wave) Model
%--- Ps: pressure at stapes
%--- Zc: input impedance of the cochlea

Vb=-P./Zp;			%--- eq. (5) and (12) in Neely's paper
Di = Hi.*Vb./s(:,ones(1,Nx));	%--- displacement of IHC (eq. 24)

%--------------------------------------------------------------------------------


%semilogx(f,db(Di(:,270:10:Nx))), title('Displacement of IHC'), xlabel('f[Hz]')


