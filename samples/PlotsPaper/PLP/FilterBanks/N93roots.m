%|------------- Neely's 1993 Model (JASA 94 (1) July 93 pp 137) -----------------------|
%|
%| Computation of the ZEROS and POLES of the partition impedance Zp(f,x)
%|
%| Input:
%|	x - row vector for x positions (must be between 0 and L=2.5cm)
%|	gama - feedback gain (default: gama=1)
%| Output:
%|	wZa;wHz;wHp - angular frequencies of quadratic factors 
%|	aZa;aHz;aHp - aZa=wZa/Qza, etc, where Q's are the quality factors
%|	wz1,wz2,wp1,wp2 - single order zeros and poles of the impedance
%|
%|------ Definition: ----------
%|
%|        Mb   (s^2+s*aZa+wZa^2)(s^2+s*aHz+wHz^2)   (s+wz1)(s+wz2) 
%|  Zp = --- * ---------------------------------- * --------------
%|        s          (s^2+s*aHp+wHp^2)              (s+wp1)(s+wp2)
%|
%|NOTE: Zp has 6 zeros and 5 poles. 2 poles and 2 zeros are always real
%|
%| function [wZa,wHz,wHp,aZa,aHz,aHp,wz1,wz2,wp1,wp2]=N93roots(x,gama)
%|
%| F. Perdigao, Coimbra, Feb. 1996 (fp@it.uc.pt)
%|===============================================================================|


function [wZa,wHz,wHp,aZa,aHz,aHp,wz1,wz2,wp1,wp2]=N93roots(x,gama)

if nargin==1, gama=1; end	%-- sets the default gama

[n,Nx]=size(x);
if (n>1),
  error('--- x must be a row ---')
end

if length(gama)>1,
 error('--- only one value for gama ---')
end

%-------------------- INITIALIZATIONS -----------------------------

rho =1;
L=2.5;
wHB=2*pi*1390;
bw=0.01;
Ap=bw*L/Nx;  %-- =bw*delta

%------------------------------------------------------------------

%---- 1st part: coeficients variation with x -------------

xx=[0,L/2,L];

Kb = exp(polyval(polyfit(xx,log([1.14e8,  4.19e6,  5.97e4 ]),2)  ,x));
Kt = exp(polyval(polyfit(xx,log([1.99e4,  2.21e4,  3.16e4 ]),2)  ,x));
K0 = exp(polyval(polyfit(xx,log([1.05e4,  9.23e3,  1.25e4 ]),2)  ,x));
Mb = exp(polyval(polyfit(xx,log([9.14e-6, 9.6e-6,  1.06e-5]),2)  ,x));
Mt = exp(polyval(polyfit(xx,log([5.64e-4, 1.02e-3, 1.06e-2]),2)  ,x));
Rb = exp(polyval(polyfit(xx,log([2.08e-2, 2.03e-2, 1.88e-2]),2)  ,x));
Rt = exp(polyval(polyfit(xx,log([149,     63.4,    27     ]),2)  ,x));
R0 = exp(polyval(polyfit(xx,log([2037,    282,     38     ]),2)  ,x));
gf = exp(polyval(polyfit(xx,log([1.42e5,  1.05e4,  3.68e2 ]),2)  ,x));
tauf=exp(polyval(polyfit(xx,log([1.4e-4,  6.92e-4, 5.29e-3]),2)  ,x));
taur=exp(polyval(polyfit(xx,log([1.35e-4, 3.61e-4, 2.5e-3 ]),2)  ,x));
%Ac= exp(polyval(polyfit(xx,log([5.52e-3, 3.17e-3, 4.27e-3]),2)  ,x));
alfa =       polyval(polyfit(xx,[0.081,   0.035,   0      ] ,2)  ,x);
gr=0.1;

%-- natural frequencies of impedances

wb = sqrt(Kb./Mb);
wt = sqrt(Kt./Mt);
w0t = sqrt(K0./Mt);
w1t = sqrt((Kt+K0)./Mt);
wc = sqrt(1 ./(tauf.*taur));

%-- coeficients of polynomials in s=j*w --

ab = Rb./Mb;
at = Rt./Mt;
a0t= R0./Mt;
a1t= (R0+Rt)./Mt;
ac = (tauf+taur)./(tauf.*taur);
An = (gf*gr)./(tauf.*taur);	%-- gain for Hc
%------------------------------------------------



for k=1:Nx,		%----- for each cochlear section... 

%-- quadratic sections

Sb = [1, ab(k), wb(k)^2];	%--- quadratic factor for Zb
St = [1, at(k), wt(k)^2];	%--- quadratic factor for Zt
S0t= [0, a0t(k),w0t(k)^2];	%--- quadratic factor for Z0
S1t= [1, a1t(k),w1t(k)^2];	%--- quadratic factor for Zt//Z0
Sc = [1, ac(k), wc(k)^2];	%--- denominator quadratic factor for Hc

%----- Zp as a ratio of polynomials in s ---------

%-- Zp = (Mb/s)*{ Sb*[S1t*Sc+gama*An*St]+(Mt/Mb)*St*Sc*p0t }/{ [S1t*Sc+gama*An*St] }

%------------------------------------------------------------------------------------

Num= conv(Mb(k)*Sb,conv(S1t,Sc)+[0,0,gama*An(k)*St]) + conv(conv(Mt(k)*St,Sc),S0t);
Den= conv([1,0],conv(S1t,Sc)+gama*An(k)*[0,0,St]);

[QSN,QSD,G(k)]=tf2qs(Num,Den);	%--- Pass to 2nd order factors. G=Mb (tested!)

aHz(k)=QSN(1,2);		%--- 2nd order factors: Sk = (s^2 + ak*s + wk2)
wHz(k)=sqrt(QSN(1,3));

aZa(k)=QSN(2,2);
wZa(k)=sqrt(QSN(2,3));

a3(k)=QSN(3,2);			%-- QSN(3,:) has only real roots
w32(k)=QSN(3,3);
r=roots([1,a3(k),w32(k)]).';	%-- These roots are always real for any x
wz1(k)=-r(1);
wz2(k)=-r(2);

aHp(k)=QSD(1,2);
wHp(k)=sqrt(QSD(1,3));

a6(k)=QSD(2,2);			%-- QSD(2,:) has only real roots
w62(k)=QSD(2,3);
r=roots([1,a6(k),w62(k)]).';	
wp1(k)=-r(1);
wp2(k)=-r(2);


end	%---- end for k=1:Nx

%-----------------------------------------------------------------------------

%---- Real roots wz2 and wp2 are (almost) always equal (min ratio=0.92 for x>2.3cm)
%---- Real roots wz1 and wp1 only affect responses for high freq. and x->L.



return

%----- TEST ------------

f=logspace(2,5,120)';
s=j*2*pi*f;
invs=1 ./s;
uf=ones(size(f));
ux=ones(size(x));
s2=(s.^2)*ux;

Zp=Neely93Z(f,x);


Z1 = (invs*Mb).*(s2+s*aZa+uf*(wZa.^2)).*(s2+s*aHz+uf*(wHz.^2))./...
                (s2+s*aHp+uf*(wHp.^2));
Z2= (s*ux+uf*wz1).*(s*ux+uf*wz2)./((s*ux+uf*wp1).*(s*ux+uf*wp2));

Zp_test=Z1.*Z2;

semilogx(f,db(Zp(:,1:Nx)),f,db(Zp_test(:,1:Nx)/Ap),':')
title('Zp vs freq.'), xlabel('f[Hz]'),ylabel('dB'), pause

plot(x,db(Ap*Zp(1:10:120,:)),'-',x,db(Zp_test(1:10:120,:)),':')
title('Zp vs distance x'), xlabel('x[cm]'),ylabel('dB'), pause

%--- wz2 and wp2
plot(x,wz2,x,wp2),title('wz2 and wp2')
