function [WZa,WHz,WHp,QZa,QHz,QHp]=KanisG(gama,x)

%+============== Kanis and de Boer Model (linear active case)  ================|
%|
%| ZKB parameter variation with gama.
%| Input:
%|	gama - vector with gama values (gama=c0/M (feedback gain) )
%| Output:
%|      WZa;WHz;WHp : angular normalized freq. for ZKB=Za*Ha
%|	QZa;QHz;QHp : quality factors for ZKB
%| Notes: ZKB scals, i.e. ZKB(s/w0) is independent of x. So, the normalized
%|        parameters do not vary either.
%|        WHp and QHp are equal to sigma and sigma/dSC.
%|
%| see KanisZ.m
%|
%| function [WZa,WHz,WHp,QZa,QHz,QHp]=KanisG(gama,x)
%|----------------------------------------------
%| F. Perdigao, Coimbra, Feb. 1996 (fp@it.uc.pt)
%+==============================================================================|


Ng=max(size(gama));

%------ Model constants ------------------
alfa = 300;
S0 = 1e10;	%-- stiffness value for x=0
M=0.5;
delta=0.4;
dSC=0.14;
sigma=0.7;
c0 = 0.06; 
%--- L=0.032;

%--------------- w0 ------------------------
wM=sqrt(S0/M);
w0=wM*exp(-alfa*x/2);
w02=w0.^2;


%------- normalized impedance with u=s/w0 --------------
%
%   zbm(u) = M*w0*(u^2 + delta*u +1)/u;
%   zohc(u)= u*(u+1)/(u^2+udSC+sigma^2)
%   zbm/(M*w0) = SBM/u
%   zohc/(M*w0)= gama*c0*SOHCn/SOHCd;
%   ZKB(u) = M*w0*( SBM/u - gama*c0*SOHCn/SOHCd ) = Za*Ha

u = [1,0];
SBM   = [1,delta,1];	%-- polynomial in u
SOHCn = [1 , 1 , 0];
SOHCd = [1, dSC, sigma^2];

WHp=sigma;
QHp=sigma/dSC;


%----- Normalized frequency and Q variation with gama: --------

for n=1:Ng,

 Num = conv(SBM,SOHCd)-gama(n)*c0*[0,conv(u,SOHCn)];

 SN=tf2qs(Num,1);

 W1=sqrt(SN(1,3)); W2=sqrt(SN(2,3));
 if (W1>W2),
  WZa(n)=W1; QZa(n)=W1/SN(1,2);
  WHz(n)=W2; QHz(n)=W2/SN(2,2);
 else
  WZa(n)=W2; QZa(n)=W2/SN(2,2);
  WHz(n)=W1; QHz(n)=W1/SN(1,2);
 end

end	%--for Ng

%--- Za = (M/s)*(s^2 + s*WZa*w0/QZa + (WZa*w0)^2)
%--- Ha = (s^2 + s*WHz*w0/QHz + (WHz*w0)^2)/((s^2 + s*WHp*w0/QHp + (WHp*w0)^2)
