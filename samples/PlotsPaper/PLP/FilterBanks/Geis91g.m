%|======================= GEISLER Model 1991  =================================|
%| Computation of ZG parameter variation with gama
%| Input:
%|	gama - column vector for required gama values
%|	x - row vector for cochlear positions
%| Output:
%|      WZa;WHz;WHp : normalized angular freq. for ZG=Za*Ha (Ng x Nx)
%|	QZa;QHz;QHp : quality factors for ZG (Ng x Nx)   
%| Notes:
%|      G ~= M*w0 has small variation with gama and x: 0.0486 < G/w0 < 0.05=M
%|      gama should not be too small (<0.1) in order to give best fittings
%|Definitions: see geislerz.m
%|
%| function [WZa,WHz,WHp,QZa,QHz,QHp]=Geis91G(gama,x)
%|----------------------------------------------
%| F. Perdigao, Coimbra, Feb. 1996 (fp@it.uc.pt)
%|=============================================================================|

function [WZa,WHz,WHp,QZa,QHz,QHp]=Geis91G(gama,x)

[n,Nx]=size(x);
if (n>1),
  error('--- x must be a row ---')
end
[Ng,n]=size(gama);
if (n>1),
  error('--- gama must be a column ---')
end

%--------- Model constants -----------------------
L=2.4;
M=0.05;
d=0.2;

%------------- w0 and alfa definition : -----------

w0=2*pi*( 456*(10 .^(2.1*(1-x/L)) - 0.8) );
alfa=(3.05./(3+exp(6*x/L-3)));

%------------ Analysis with normalized frequency u=jw/w0 -------------

wn=logspace(-2,1,200)';		%-- Normalized frequency in 3 decades
u=j*wn;
u2=-(wn.^2);
wt=ones(200,1);
n=fix(3*200/4)+1;
wt(151:200)=zeros(1,50);	%--- wt: freq. values for which z are not modeled
%%kass = fix(Nx/10);


%---- ZG=ZBM+ZOHC; ZOHC=(Kp/s)/(1+F/Ks)
%---- ZBM = s*M+R+Kb/s; Kb=M*w0^2; R=d*Kb/w0; d=0.2; ZBM=(M/s)(s^2+s*d*w0+w0^2)
%---- Kp=Kr*Kc/Ks; Ks=Kr+Kc; Kr=Kc=0.9Kb; Kp=0.45*Kb=0.45*M*w0^2
%---- F=-j*gama*Kb*alfa*exp(-s*T); F/Ks= -j*(gama/1.8)*alfa*exp(-s*T); T=pi/w0
%---- results in (with s=j*w; u=s/w0);
%---- ZG(s) = (M/s)*(s^2+s*d*w0+w0^2 + 0.45*w0^2/(1-j*(gama/1.8)*alfa*exp(-s*pi/w0)))
%---- ZG(u) = (M*w0/u)*(u^2+u*d+1 + 0.45/(1-j*(gama/1.8)*alfa*exp(-u*pi)))

zbm = ((u2+u*d+1)./u) * (M*w0);

for n=1:Ng,

 Ff = (-j*(gama(n)/1.8)*exp(-u*pi)) * alfa;
 zohc= ( (1 ./u) * (0.45*M*w0) )./(1+Ff);
 z=zbm+zohc;

%------------- Modelization of z as: ZGmod = Za*Ha  ------------------

for k=1:Nx,

[B,A]=invfreqs(z(:,k),wn,4,3,wt);	%-- Rational Z with 4 zeros and 3 poles

%---- assessment of fitting (already tested)---------

%if ~rem(k-1,kass),
%hz=freqs(B,A,wn);
%semilogx(wn,db(hz),wn,db(z(:,k))), title(['n=',int2str(n),' k=',int2str(k)]), pause
%semilogx(wn,fase(hz),wn,fase(z(:,k))), pause
%end

[SN,SD,G(n,k)]=tf2qs(B,A);	%--- pass to quadratic sections  (G=M*w0, almost)

%--------- parameters ---------------
a1(n,k)=SN(1,2);
W12(n,k)=SN(1,3);
a2(n,k)=SN(2,2);
W22(n,k)=SN(2,3);
a3(n,k)=SD(1,2);
W32(n,k)=SD(1,3);
%---WZp(n,k)=SD(2,3);	%--- This freq. is always smaller than 3e-3,
			%--- so, we can take it as zero (a pole at origin)
SD(2,3)=0; 		%--- a pole is inserted at origin

%--------- tests ---------------
%Num=G(n,k)*conv(SN(1,:),SN(2,:));	%-- numerator (zero freq. pole inserted)
%Den=conv([1,0],SD(1,:));		%-- denominator (zero freq. pole inserted)
%Zp=freqs(Num,Den,wn);
%semilogx(wn,db(Zp),wn,db(z(:,k))), title(['model vs. z gama=',num2str(gama(n))]), pause

end
end

%-------------- Non-normalized parameters -----------

WZa=sqrt(W12);
QZa=sqrt(W12)./a1;
WHz=sqrt(W22);
QHz=sqrt(W22)./a2;
WHp=sqrt(W32);
QHp=sqrt(W32)./a3;

