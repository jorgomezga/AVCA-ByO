  function [P,Zin,Zeq]=TlinePs(Ps,Zs,Zp)
% function [P,Zin,Zeq]=TlinePs(Ps,Zs,Zp)
% Computs cochlear pressures and the input impedance of the 
% equivalent transmission line cochlear circuit (long wave approx.).
% Input: Pressure at stapes
%
% Input: Partition impedances Zp(Nf,Nx)	, Nf,Nx -> no. of points in f,x;
%	 Series impedances Zs: Zs(Nf,1) or Zs(Nf,Nx)
%	 Ps - Fluid pressure at stapes
% Output: Pressures P(1)..P(N)
%	 Input impedances Zin=Zeq(1); Zeq
%------------------------------------------------------------------------------
%
% Example with N=5:
%
%         Zs1  P1  Zs2  P2  Zs3  P3  Zs4  P4  Zs5  P5
%     --@@@@--|--@@@@--|--@@@@--|--@@@@--|--@@@@--+
%    +        #        #        #        #        #    
%   Ps    Zp1 #    Zp2 #    Zp3 #    Zp4 #    Zp5 #
%    -        # |->    #        #        #        #
%     ________|_|______|________|________|________+ 
%   x=0     x=Dx|     x=2Dx    x=3Dx    x=4Dx    x=5Dx=Lx
%		|
%		Zeq(2)
% Ps=P(x=0); P6=P(N)=P(x=L)
%
% Related funcions TLinePs2 TLineVs TLineME TlineME2 
%------------------------------------------------------------------------------

% Fernando Perdigao, Coimbra, 1996 (fp@dee.uc.pt)

[Nf,Nx]=size(Zp);
[N1,N2]=size(Zs);

if (Nf~=N1), error('Size error in Zs'), end
if N2==1, Zs=Zs(:,ones(1,Nx)); end		%-- if Zs is constant in x
[N1,N2]=size(Ps);
if length(Ps)>1,
 if ((N1~=Nf)|(N2>1)), error('Size error in Ps'), end;
 Ps=Ps(:,ones(1,Nx));
end

%-- Equivalent impedances ---
%----------------------------

Zeq(:,Nx)=Zs(:,Nx)+Zp(:,Nx);	%-- at helicotrema

for k=Nx-1:-1:1,
 Zeq(:,k) = Zs(:,k) + Zp(:,k).*Zeq(:,k+1) ./ ( Zp(:,k)+Zeq(:,k+1) );
end

Zin=Zeq(:,1);		%-- coclhear input impedance (at x=0)

%--- Pressures computation with equivalent impedances
%----------------------------------------------------

P=Ps.*(cumprod((1-Zs./Zeq).').');		%-- P(k)=P(k-1)*(1-Zs(k)/Zeq(k))

%--------------------------------------------------------------