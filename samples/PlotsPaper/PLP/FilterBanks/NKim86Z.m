function [Zp, Zpmod, fQ, Mef] = NKim86Z(f, x, gamma)
% function [Zp, Zpmod, fQ, Mef] = NKim86Z(f, x, gamma)
% Neely and Kim's, 1986 model (JASA 79(5) May 86)
% Computation of partition impedance Zp(f,x)
% Input:
% f - column vector for required frequencies
% x - row vector for x positions
% gamma - feedback gain (default: gamma=1)
% Output:
% Zp - Model partition impedance((Nf x Nx) elements)
% Zpmod - Partition impedance as: Zpmod=Za*Ha (Nf x Nx)
% fQ - Parameters of Zpmod: frequencies and quality factors:
% fQ=[fZa;fHz;fHp;QZa;QHz;QHp] (matrix of 6 x Nx elements)
%
%                Mef                         (s^2+s*wHz/QHz+wHz^2)
%       Zpmod =  ---*(s^2+s*wZa/QZa+wZa^2) * ---------------------
%                 s                          (s^2+s*wHp/QHp+wHp^2)
%
%                |<---------Za----------->| |<---------Ha-------->|
%
%	Mef = (g/b)*m1 =0.0075 : equivalent mass

%----------------------------------------------
% F. Perdigao, Coimbra, Feb. 1996 (fp@it.uc.pt)
%----------------------------------------------


if nargin==2, gama=1; end

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

%-----------------------------------------------------------------
%-------------------- INITIALIZATIONS -----------------------------

%h=0.1;
%km=2.1e5;
%cm=400;
%mm=45e-3;
%As=0.01;
%Am=0.35;
%gm=0.5;
%rho=1;
%L=2.5;
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


s=j*2*pi*f;
invs=1 ./s;
s2 = -(2*pi*f).^2;
umf=ones(Nf,1);
umx=ones(1,Nx);

Z1=invs*k1 + c1(umf,:) + m1*s(:,umx);   %--- Z1=k1/s+c1+s*m1;
Z2=invs*k2 + c2(umf,:) + m2*s(:,umx);	%--- Z2=k2/s+c2+s*m2;
Z3=invs*k3 + c3(umf,:);			%--- Z3=k3/s+c3;
Z4=invs*k4 + c4(umf,:);			%--- Z4=k4/s+c4;
Hc=Z2./(Z2+Z3);

Zp=(g/b)*( Z1 + Hc.*(Z3-gama*Z4) );

if nargout==1, return, end

%------Caracterization of Zp as:  Zpmod=Za*Ha  -------------


%-- Z1 = (m1/s)*(s^2 + s*a1*w1 + w1^2) ;   a1=c1/m1;    w1^2 = k1/m1;
%   Z2 = (m2/s)*(s^2 + s*a2*w2 + w2^2) ;   a2=c2/m2;    w2^2 = k2/m2;
%   Z3 = c3 + k3/s;    Z4=c4 + k4/s;
%   Hc = Z2/(Z2+Z3) = (s^2 + s*a2 + w2^2)/(s^2 + s*a5 + w5^2); 
%   a5 = (c2+c3)/m2;   w5^2=(k2+k3)/m2;
%   Zp = (g/b)*(m1/s)*( s^2 + s*a1 + w1^2  + Hc*(s*a6 + w6^2) ) = Za*Ha
%   a6 = (c3-gama*c4)/m1  ;  w6^2 = (k3-gama*k4)/m1;

Mef=(g/b)*m1;
a1=c1/m1; 
a2=c2/m2;
a5=(c2+c3)/m2;
a6=(c3-gama*c4)/m1;
w12=k1/m1;
w22=k2/m2;
w52=(k2+k3)/m2;
w62=(k3-gama*k4)/m1;


for k=1:Nx,

SBM = [1,a1(k),w12(k)];
S2  = [1,a2(k),w22(k)];
S5  = [1,a5(k),w52(k)];
S6  = [0,a6(k),w62(k)];

Num = Mef*(conv(SBM,S5)+conv(S2,S6));


SN=tf2qs(Num,1);	%--- pass to 2 quadratic sections

%--------- parameters ---------------

wZa(k) = sqrt(SN(1,3));
QZa(k) = sqrt(SN(1,3))/SN(1,2);
wHz(k) = sqrt(SN(2,3));
QHz(k) = sqrt(SN(2,3))/SN(2,2);


end

wHp = sqrt(w52);
QHp = sqrt(w52)./a5;

Za= Mef*( s(:,umx) +umf*(wZa./QZa) + invs*(wZa.^2) );
Ha=     (s2(:,umx) +  s*(wHz./QHz) +  umf*(wHz.^2) )./ ...
        (s2(:,umx) +  s*(wHp./QHp) +  umf*(wHp.^2) );
Zpmod = Za.*Ha;

%--------- test ---------------
%semilogx(f,db(Zp(:,50)),f,db(Zpmod(:,50)))


%----- curve crossing -----------

dQ=abs(QZa(2:Nx)-QZa(1:Nx-1));
i=find(dQ>1);	%--- too much change in QZa?

%---- exchange zeros and quality factors between Za and Ha:
if length(i)==1,
   i=i+1; 
   w=wZa(i:Nx); wZa(i:Nx)=wHz(i:Nx); wHz(i:Nx)=w;
   Q=QZa(i:Nx); QZa(i:Nx)=QHz(i:Nx); QHz(i:Nx)=Q;
end


%--- The 6 parameters in matrix form (7xNx)
tp=2*pi;
fQ = [wZa/tp;wHz/tp;wHp/tp;QZa;QHz;QHp]; 

