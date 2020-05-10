function [ZG,ZGmod,WQ,M]=Geis93Z(f,x,gama)
%------------- GEISLER's Model (1993)  -----------------------|
%
% Computation of partition impedance ZG(f,x)
% and related simplified impedance: ZGmod = Za*Ha
% Input:
%	f - column vector of response's desired frequencies
%	x - row vector of cochlear positions
%	gama - feedback gain (default: gama=1.29)
% Output:
%	ZG - Partition Impedance of Geisler's model ((Nf x Nx) elements)
%	ZGmod - Alternative impedance ZGmod=Za*Ha (Nf x Nx)
%	WQ  - normalized angular frequencies and quality factors of ZGmod
%		WQ=[WZa,WHz,WHp,QZa,QHz,QHp] (6 elements)
%    Definitions: 
%	ZG=ZBM+ZOHC; ZOHC=(Kp/s)/(1+F/Ks)
%	ZBM = s*M+R+Kb/s; Kb=M*w0^2; R=d*Kb/w0; ZBM=(M/s)(s^2+s*d*w0+w0^2)
%	Kp=Kr*Kc/Ks; Ks=Kr+Kc; Kr=1.17Kb; Kc=0.59Kb; 
%	Kp=0.392216*Kb = rpb*Kb; Ks=1.76*Kb=rsb*Kb
%	Ffb= gama*Kb*[(wr-s)/(wr+s)]*exp(-s*T); T=n*pi/w0
%
%	With u=jw/w0, ZG/w0 scals, i.e. is independent of x.
%
% Reference: Hearing Research, 68(1993) 253-262
%
%	Alternate impedance: ZGmod=Za*Ha
%
%                 M                            (s^2+s*wHz/QHz+wHz^2)
%       ZGmod =  --- * (s^2+s*wZa/QZa+wZa^2) * ---------------------
%                 s                            (s^2+s*wHp/QHp+wHp^2)
%
%               |<-----------Za------------>| |<---------Ha--------->|
%
% function [ZG,ZGmod,WQ,M]=Geis93Z(f,x,gama)
%----------------------------------------------
% F. Perdigao, Coimbra, Feb. 1996 (fp@it.uc.pt)
%------------------------------------------------------------------------------|

if nargin==2, gama=1.29; end	%--- default gama, near the optimum value

[n,Nx]=size(x);
if (n>1),
  error('--- x must be a row ---')
end
[Nf,n]=size(f);
if (n>1),
  error('--- f must be a column ---')
end

if length(gama)>1,
 error('--- only one value for gama ---')
end


%--------- Model Constants and definitions ---------------------
M=0.05;		%--- gr/cm^2
L=2.4; 		%--- cm
n=0.4;		%--  T=n*pi/w0
d=0.6;		%--  d=R/(M*w0) = R*w0/Kb

%-- D = C + cos(psi)		(eq. 4)
%-- C=(p/b)*cotan(csi)		(eq. A9, A10)
%-- DA/B =(C+cos(psi))*(A/B)	A/B=1.71
%-- p/b=0.57 ; A/B=1.71; psi=18o; csi=32o
%-- C= 0.57/tan(32*pi/180);

DAB = 1.71*( 0.57/tan(32*pi/180) + cos(18*pi/180) );	%--- = DA/B

%--- rpb = Kp/Kb = (Kr//Kc)/Kb = 1.17*0.59/(1.17+0.59)
%--- rsb = Ks/Kb = (Kr+Kc)/Kb = 1.76
rpb = 1.17*0.59/(1.17+0.59);
rsb = 1.76;
%------------------------------------------------------------

%------ Cochlear Map: f0: 57KHz downto 90Hz -----------------

w0=2*pi*( 456*(10 .^(2.1*(1-x/L)) - 0.8) );
Kb = M*(w0.^2);
Kr=1.17*Kb; 
Kc=0.59*Kb; 
R=d*Kb./w0;

T=n*pi./w0;
F0=gama*Kb;				
wr=w0/20;

%------------ PARTITION IMPEDANCE Z ----------------------------------

s=j*2*pi*f;
invs = 1 ./s;
uf=ones(size(f));
ux=ones(size(x));


ZBM = s*(M*ux) + uf*R + invs*Kb;
apf = (invs*wr -1)./(invs*wr+1);	%--- All-Pass Filter
Ffb = (uf*F0).*apf.*exp(-s*T);

ZOHC = DAB*(invs*(Kr.*Kc))./(uf*(Kr+Kc) + Ffb);

ZG = ZBM + ZOHC;

if nargout==1, return, end

%=====================================================================

%------------ Analysis with normalized frequency u=jw/w0 -------------
%---- ZG=ZBM+ZOHC; ZOHC=(Kp/s)/(1+F/Ks)
%---- ZBM = s*M+R+Kb/s; Kb=M*w0^2; R=d*Kb/w0; ZBM=(M/s)(s^2+s*d*w0+w0^2)
%---- Kp=Kr*Kc/Ks; Ks=Kr+Kc; Kr=1.17Kb; Kc=0.59Kb; 
%---- Kp=0.392216*Kb = rpb*Kb; Ks=1.76*Kb=rsb*Kb
%---- Ffb= gama*Kb*[(wr-s)/(wr+s)]*exp(-s*T); T=n*pi/w0

%---- results in: (with s=j*w; u=s/w0)
%---- ZBM(u) = (M*w0/u)*(u^2+u*d+1)
%---- F0=gama*Kb = gama*M*w0^2; wr=w0/20 => (wr-s)/(wr+s) = (1/20-u)/(1/20+u);
%---- ZOHC(u)= (M*w0/u)*(DAB*rpb)/( 1 + (gama/rsb)*(0.05-u)/(0.05+u)*exp(-u*n*pi) )

%---- THEN:

%---- z = ZG/(M*w0/u) is independent of x !!! ---------------------------------

Nu=500;			%--- Number of points in u or wn
wn=logspace(-2,0.2,Nu)';
u=j*wn;
u2=-wn.*wn;
uw=ones(size(wn));

ffbKs = (gama/rsb)*((0.05-u)./(0.05+u)).*exp(-u*n*pi);	%--- Ffb/(Kr+Kc)
z = u.^2+u*d+1 + (DAB*rpb)./(1+ffbKs);

AH=[u.*u2,u2,u,uw,-u.*z,-z];		%----- LMS fitting with 4 zeros and 2 poles.
cH=-u2.*u2 + z.*u2;
A=real(AH'*AH); cH=real(AH'*cH);
cH=A\cH;
B=[1,cH(1:4)'];
A=[1,cH(5:6)'];

[SN,SD,G]=tf2qs(B,A);

WZa=sqrt(SN(1,3));		%--- WZa2 = 1.7288
QZa=WZa/SN(1,2);		%--- aZa  = 0.1374
WHz=sqrt(SN(2,3));		%--- WHz2 = 0.0791
QHz=WHz/SN(2,2);		%--- aHz  = 0.8489
WHp=sqrt(SD(1,3));		%--- WHp2 = 0.0793
QHp=WHp/SD(1,2);		%--- aHp  = 0.2491
%----WZa2*WHp2/WHp2;		%--- 1.7288 

%---- test ---------------------------
%wn=logspace(-3,2,Nu)'; u=j*wn;
%ffbKs = (gama/rsb)*((0.05-u)./(0.05+u)).*exp(-u*n*pi);	%--- Ffb/(Kr+Kc)
%z = u.^2+u*d+1 + (DAB*rpb)./(1+ffbKs);
%hz=freqs(B,A,wn); 
%semilogx(wn,db(hz),wn,db(z)), pause, semilogx(wn,fase(hz),wn,fase(z))
%plot(1:Nu,db(z),1:Nu,db(hz),1:Nu,db(z-hz)), pause


%---- fZa = 599.275*(10^(2.1*(1-x/L)) - 0.8)	= sqrt(W12)*w0/(2*pi)
%---- fHz = 128.608*(10^(2.1*(1-x/L)) - 0.8)	= sqrt(W22)*w0/(2*pi)
%---- fHp = 128.513*(10^(2.1*(1-x/L)) - 0.8)	= sqrt(W32)*w0/(2*pi)

WQ = [WZa,WHz,WHp,QZa,QHz,QHp];

Za = M*( s*ux + uf*(WZa*w0/QZa) + invs*((WZa*w0).^2) );
Ha = ( (s.^2)*ux + s*(WHz*w0/QHz) + uf*((WHz*w0).^2) )./ ...
     ( (s.^2)*ux + s*(WHp*w0/QHp) + uf*((WHp*w0).^2) );
ZGmod = Za.*Ha;

