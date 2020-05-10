  function ZL=N93ZL(f,gama,Ap)
% function ZL=N93ZL(f,gama,Ap);
%
%   N93ZL: Computs helicotrema impedance for no reflections in Neely's model
%           (same case as if an infinite cochlea)
%   Input: f - freq. column vector
%	   gama- feedback coeficient
%   Output: ZL - "helicotrema" impedance  



%---- consts: -----------
Nf=length(f);
if nargin==1, gama = 1; end
L=2.5;
rho=1;
Dx=L/500;
if nargin <3, bw=0.01; Ap=bw*Dx; end


N=300;
Lmax=4;
x1=linspace(L,Lmax,N);%-- Dx is the same: 2.5/500=3/600=3.5/700=4/800
s=j*2*pi*f;


%Ac=exp(polyval(polyfit([0,L/2,L],log([5.52e-3, 3.17e-3, 4.27e-3]),2)  ,2.5))
c=polyfit([0,L/2,L],log([5.52e-3, 3.17e-3, 4.27e-3]),2);
b=[L^2,L,1;2*L,1,0;Lmax^2,Lmax,1]\[polyval(c,L);2*c(1)*L+c(2);polyval(c,L/2)];
Ac=exp(polyval(b,x1));


Mf=2*rho*Dx./Ac;
Rf=2e2;
zs = Rf+s*Mf;		%--- series impedance
zp=Neely93Z(f,x1)/Ap;	%--- parallel impedance
zp(:,N)=sqrt(zp(:,N).*zs(:,N));


[ZinN,Zr]=coc_zin(zs,zp);	%--- gets the "left" impedance at every point

ZL=Zr(:,1)-zs(:,1);		%-- sets the "helicotrema impedance at x=L=2.5 cm

%semilogx(f,db(Zr(:,1:20:N)))