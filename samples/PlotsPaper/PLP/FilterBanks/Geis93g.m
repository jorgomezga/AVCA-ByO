%|======================= GEISLER-1993 Model =================================|
%| Computation of ZG's parameter variation with gama
%| Input:
%|	gama - column vector for required gama values
%|
%| Output:
%|      WZa;WHz;WHp : normalized angular freq. for ZG=Za*Ha (Ng x 1)
%|	QZa;QHz;QHp : quality factors for ZG (Ng x 1)   (ZG scals)
%|Definitions: see Geis93Z.m
%|
%| function [WZa,WHz,WHp,QZa,QHz,QHp]=Geis93G(gama)
%|----------------------------------------------
%| F. Perdigao, Coimbra, Feb. 1996 (fp@it.uc.pt)
%|=============================================================================|

function [WZa,WHz,WHp,QZa,QHz,QHp]=Geis93G(gama)

Ng=length(gama);

%--------- Model Constants ---------------------
M=0.05;		%-- gr/cm^2
n=0.4;		%-- T=n*pi/w0
d=0.6;		%-- d=R/(M*w0) = R*w0/Kb
		%-- D = C + cos(psi)		(eq. 4)
		%-- C=(p/b)*cotan(csi)		(eq. A9, A10)
		%-- DA/B =(C+cos(psi))*(A/B)	A/B=1.71
		%-- p/b=0.57 ; A/B=1.71; psi=18o; csi=32o
		%-- C= 0.57/tan(32*pi/180);
		%-- rpb = Kp/Kb = (Kr//Kc)/Kb = 1.17*0.59/(1.17+0.59)
		%-- rsb = Ks/Kb = (Kr+Kc)/Kb = 1.76
DAB = 1.71*( 0.57/tan(32*pi/180) + cos(18*pi/180) );	%--- = DA/B
rpb = 1.17*0.59/(1.17+0.59);
rsb = 1.76;
%--------------------------

Nu=500;
wn=logspace(-3,0.25,Nu)';
u=j*wn;
u2=-wn.*wn;
u3=u.*u2;
u4=u2.*u2;
uw=ones(size(wn));


ffbKs = (1/rsb)*((0.05-u)./(0.05+u)).*exp(-u*n*pi);	%--- Ffb/(Kr+Kc) with gama=1


for k=1:Ng,	%--- see Geis93Z.m

z = u2+u*d+1 + (DAB*rpb)./(1+gama(k)*ffbKs);
AH=[u3,u2,u,uw,-u.*z,-z];		%----- LMS fitting with 4 zeros and 2 poles.
cH=-u4+z.*u2;
A=real(AH'*AH); cH=real(AH'*cH);
cH=A\cH;
B=[1,cH(1:4)']; A=[1,cH(5:6)'];

[SN,SD,G]=tf2qs(B,A);


WZa(k)=sqrt(SN(1,3));
WHz(k)=sqrt(SN(2,3));
WHp(k)=sqrt(SD(1,3));

QZa(k)=WZa(k)/SN(1,2);
QHz(k)=WHz(k)/SN(2,2);
QHp(k)=WHp(k)/SD(1,2);

end


