%|------------- GEISLER's Model 1991  -----------------------|
%|
%| Computation of partition impedance ZG(f,x)
%| and related simplified impedance: ZGmod = Za*Ha
%| Input:
%|	f - column vector of response's required frequencies
%|	x - row vector of cochlear positions
%|	gama - feedback gain (default: gama=1.4)
%| Output:
%|	ZG - Partition Impedance of Geisler's model ((Nf x Nx) elements)
%|	ZGmod - Alternative impedance ZGmod=Za*Ha (Nf x Nx)
%|	fQ - Parameters of ZGmod: frequencies and quality factors:
%|		fQ=[fZa;fHz;fHp;QZa;QHz;QHp] (matrix of 6 x Nx elements)
%|
%|------ Definitions: ----------
%|		ZG=ZBM+ZOHC;   ZOHC=(Kp/s)/(1+F/Ks)
%|		ZBM = s*M + R + Kb/s;
%|		Kb=M*w0^2; R=d*Kb/w0; d=0.2 => ZBM=(M/s)(s^2+s*d*w0+w0^2)
%|		Kp=Kr*Kc/(Kr+Kc); Kr=Kc=0.9Kb; Kp=0.45*Kb=0.45*M*w0^2
%|		alfa = 3.05/(3+exp(6*x/L-3))   L=2.4cm
%|		F0 = gama*Kb*alfa
%|		F=-j*F0*exp(-s*T); F/Ks= -j*(gama/1.8)*alfa*exp(-s*T); T=pi/w0
%|
%|                 M                            (s^2+s*wHz/QHz+wHz^2)
%|       ZGmod =  --- * (s^2+s*wZa/QZa+wZa^2) * ---------------------
%|                 s                            (s^2+s*wHp/QHp+wHp^2)
%|
%|               |<-----------Za------------>| |<---------Ha--------->|
%|
%| function [ZG,ZGmod,fQ,M]=Geis91Z(f,x,gama)
%|----------------------------------------------
%| F. Perdigao, Coimbra, Feb. 1996 (fp@it.uc.pt)
%+------------------------------------------------------------------------------|


function [ZG,ZGmod,fQ,M]=Geis91Z(f,x,gama)

if nargin==2, gama=1.4; end

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


%--------- Model constants -----------------------
L=2.4;		%-- CGS Units
M=0.05;
d=0.2;

%------------- w0 and alfa definition : -----------

w0=2*pi*( 456*(10 .^(2.1*(1-x/L)) - 0.8) );
alfa=(3.05./(3+exp(6*x/L-3)));


%------------- The non-normalized impedance -----------------

s=j*2*pi*f; 
s2=-(2*pi*f).^2;

zbm = M*s(:,ones(1,Nx)) + d*M*w0(ones(1,Nf),:) + (M./s)*(w0.^2);
Ff = ((-j*gama/1.8)*exp((-pi*s)*(1 ./w0))) .* alfa(ones(1,Nf),:);
zohc= ( (1 ./s)*(0.45*M*(w0.^2)) )./(1+Ff);

ZG=zbm+zohc;

if nargout==1, return, end

%------------ Analysis with normalized frequency u=jw/w0 -------------


wn=logspace(-2,1,200)';		%-- Normalized frequency in 3 decades and 200 points
u=j*wn;
u2=-(wn.^2);

%---- ZG=ZBM+ZOHC; ZOHC=(Kp/s)/(1+F/Ks)
%---- ZBM = s*M+R+Kb/s; Kb=M*w0^2; R=d*Kb/w0; d=0.2; ZBM=(M/s)(s^2+s*d*w0+w0^2)
%---- Kp=Kr*Kc/Ks; Ks=Kr+Kc; Kr=Kc=0.9Kb; Kp=0.45*Kb=0.45*M*w0^2
%---- F=-j*gama*Kb*alfa*exp(-s*T); F/Ks= -j*(gama/1.8)*alfa*exp(-s*T); T=pi/w0
%---- results in (with s=j*w; u=s/w0);
%---- ZG(s) = (M/s)*(s^2+s*d*w0+w0^2 + 0.45*w0^2/(1-j*(gama/1.8)*alfa*exp(-s*pi/w0)))
%---- ZG(u) = (M*w0/u)*(u^2+u*d+1 + 0.45/(1-j*(gama/1.8)*alfa*exp(-u*pi)))

 zbm = ((u2+u*d+1)./u) * (M*w0);		%--- matrix Nf x Nx
 Ff = (-j*(gama/1.8)*exp(-u*pi)) * alfa;	%--- matrix Nf x Nx
 zohc= ( (1 ./u) * (0.45*M*w0) )./(1+Ff);	%--- matrix Nf x Nx
 z=zbm+zohc;					%--- matrix Nf x Nx


%------------- Modelization of z as: ZGmod = Za*Ha  ------------------

wt=ones(200,1);
%--n=fix(3*Nf/4)+1;
wt(151:200)=zeros(1,50);	%--- wt: freq. values for which z are not modeled

%kass = fix(Nx/10);

for k=1:Nx,

[B,A]=invfreqs(z(:,k),wn,4,3,wt);	%-- Rational Z with 4 zeros and 3 poles

[SN,SD,G(k)]=tf2qs(B,A);	%--- pass to quadratic sections  (G=M*w0)
SD(2,3)=0; 			%--- a pole is inserted at origin

%---- assessment of fitting for 10 points ---------
%if ~rem(k-1,kass),
%hz=freqs(B,A,wn);
%wn1=sqrt(SN(1,3)); wn2=sqrt(SN(2,3)); wn3=sqrt(SD(1,3));
%semilogx(wn,db(hz),wn,db(z(:,k)),[wn1,wn1],[40,140],[wn2,wn2],[40,140],[wn3,wn3],[40,140]), pause
%semilogx([wn1,wn1],[40,140],[wn2,wn2],[40,140],[wn3,wn3],[40,140]), pause
%semilogx(wn,fase(hz),wn,fase(z(:,k))), pause
%end



%--------- parameters ---------------
a1(k)=SN(1,2);
W12(k)=SN(1,3);
a2(k)=SN(2,2);
W22(k)=SN(2,3);
a3(k)=SD(1,2);
W32(k)=SD(1,3);


%--------- tests ---------------
%Num=G(k)*conv(SN(1,:),SN(2,:));		%-- numerator
%Den=conv([1,0],SD(1,:));			%-- denominator
%Zp=freqs(Num,Den,wn);
%semilogx(wn,db(Zp),wn,db(z(:,k))), title('model vs. z'), pause		%-- ok

end

%-------------- Non-normalized parameters -----------

wZa=sqrt(W22).*w0;
QZa=sqrt(W22)./a2;

wHz=sqrt(W12).*w0;
QHz=sqrt(W12)./a1;

wHp=sqrt(W32).*w0;
QHp=sqrt(W32)./a3;

um=ones(Nf,1);

Za= M*( s(:,ones(1,Nx)) +um*(wZa./QZa) + (1 ./s)*(wZa.^2) );
Ha=  ( s2(:,ones(1,Nx)) + s*(wHz./QHz) + um*(wHz.^2) )./ ...
     ( s2(:,ones(1,Nx)) + s*(wHp./QHp) + um*(wHp.^2) );

ZGmod = Za.*Ha;

%--- The 6 parameters in matrix form (7xNx)
tp=2*pi;
fQ = [wZa/tp;wHz/tp;wHp/tp;QZa;QHz;QHp]; 
