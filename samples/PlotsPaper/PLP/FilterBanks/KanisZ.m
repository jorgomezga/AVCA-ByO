function [ZKB,WQ,M]=KanisZ(f,x,gama)

%|------------- Kanis and de Boer Model (linear active case)  -----------|
%|
%| Computation of partition impedance ZKB(f,x)
%| Input:
%|	f - column vector of response's desired frequencies
%|	x - row vector of cochlear positions
%|	gama - feedback gain
%| Output:
%|	ZKB - Partition Impedance of Kanis-Boer model ((Nf x Nx) elements)
%|	WQ  - normalized angular frequencies and quality factors of ZKB
%|		WQ=[WZa,WHz,WHp,QZa,QHz,QHp] (6 elements)
%|	Note:	ZKB scals, i.e. ZKB(s/w0) is independent of x
%|		ZKB angular freq.s are: WZa*w0, etc.
%|------ Definitions: ----------
%|		ZBM = (M/s)*(s^2 + s*delta*w0 + w0.^2);
%|		ZOHC = gama*c0*(s*w0*(s+w0)/(s^2 + s*dSC*w0 + (sigma*w0)^2);
%|		ZKB = ZBM-ZOHC   ;  gama*M=c0
%|		w0 = wM*exp(-alfa*x); wM=sqrt(K(x=0)/M)
%|		ZKB = Za*Ha;
%|          M
%|   Za =  --- * (s^2 + s*WZa*w0/QZa + (WZa*w0)^2)
%|          s
%|
%|        (s^2 + s*WHz*w0/QHz + (WHz*w0)^2)
%|   Ha = ----------------------------------
%|        (s^2 + s*WHp*w0/QHp + (WHp*w0)^2)
%|
%| function [ZKB,WQ,M]=KanisZ(f,x,gama)
%|
%| References: JASA 94(6), pp. 3199-3206; JASA 96(4), pp. 2156-2165
%|
%|----------------------------------------------
%| F. Perdigao, Coimbra, Feb. 1996 (fp@it.uc.pt)
%+------------------------------------------------------------------------------|



if nargin==2, gama=1; end	%-- Takes maximum gama if not specified

[n,Nx]=size(x);
if (n>1),
  error('--- x must be a row ---')
end
[Nf,n]=size(f);
if (n>1),
  error('--- f must be a column ---')
end

if length(gama)>1,
 error('--- must be only one value of gama ---')
end

%------ Model constants ------------------
alfa = 300;	%-- MKS units
S0 = 1e10;	%-- stiffness value for x=0
M=0.5;
delta=0.4;
dSC=0.14;
sigma=0.7;
c0=0.06;
%--- L=0.032;	<--------- cochlear length in [m]

%--------------- w0 ------------------------
wM=sqrt(S0/M);			%-- Maximum freq. (at windows)
w0=wM*exp(-alfa*x/2);		%-- ressonance frequencies
w02=w0.^2;

w=2*pi*f;
s=j*w; s2=-w.^2;

ux=ones(1,Nx); 
uf=ones(1,Nf);

ZBM = M*s(:,ux) + delta*M*w0(uf,:) + (M./s)*(w0.^2);

%---- c0=gama*M -------

ZOHC = (gama*c0)*( s2*w0 + s*w02 )./( s2(:,ux) + dSC*s*w0 + (sigma^2)*w02(uf,:) );
ZKB = (ZBM - ZOHC);


if nargout==1, return, end	

%------- normalized impedances with u=s/w0 --------------
%
%   zbm(u) = M*w0*(u^2 + delta*u +1)/u;
%   zohc(u)= u*(u+1)/(u^2+udSC+sigma^2)
%   zbm/(M*w0) = SBM/u
%   zohc/(M*w0)= gama*SOHCn/SOHCd;
%   ZKB(u) = M*w0*( SBM/u - gama*SOHCn/SOHCd ) = Za*Ha

u = [1,0];
SBM   = [1,delta,1];	%-- polynomial in u
SOHCn = [1 , 1 , 0];
SOHCd = [1, dSC, sigma^2];

Num = conv(SBM,SOHCd)-gama*c0*[0,conv(u,SOHCn)];
Den = conv(u,SOHCd);

%----- Normalized frequencies: --------

[SN,SD,G]=tf2qs(Num,Den);	%--- passes to quadratic factors


W1=sqrt(SN(1,3)); 
W2=sqrt(SN(2,3));
if (W1>W2),			%-- sorts frequencies: WZa>WHz
 WZa=W1; QZa=W1/SN(1,2);
 WHz=W2; QHz=W2/SN(2,2);
else
 WZa=W2; QZa=W2/SN(2,2);
 WHz=W1; QHz=W1/SN(1,2);
end

WHp=sqrt(SD(1,3));		%-- WHp = sigma = 0.7
QHp=sqrt(SD(1,3))/SD(1,2);	%-- QHp = sigma/dSC = 5

WQ = [WZa,WHz,WHp,QZa,QHz,QHp];


%------- Za and Ha in ZNK=Za*Ha, may be computed as follows:---------------
%--- Za = (M/s)*(s^2 + s*WZa*w0/QZa + (WZa*w0)^2)
%--- Ha = (s^2 + s*WHz*w0/QHz + (WHz*w0)^2)/((s^2 + s*WHp*w0/QHp + (WHp*w0)^2)


