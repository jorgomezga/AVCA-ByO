  function [P,Ps,Zin,Zeq]=TlineME(Pin,Zm,Zs,Zp)
% function [P,Ps,Zin,Zeq]=TlineME(Pin,Zm,Zs,Zp)
%
% Computes cochlear pressures and the input impedance of the 
% equivalent transmission line cochlear circuit (long wave approx.),
% with middle-ear connected.
%
% Input: series and parallel distributed impedances:
%	  series impedance: Zs: Zs(Nf,1) or Zs(Nf,Nx), Nf: no. of points in f.
%	  parallel impedances: Zp(Nf,Nx),              Nx: no. of points in x.
%	 Zm:	Middle ear impedance referred to stapes velocity
%	 Pin:	Input pressure/force    "	"	"
% Output: Pressures P(1)..P(N); Ps - pressure at stapes
%	  Input impedances Zin=Zeq(1); Zeq
%------------------------------------------------------------------------------
%
% Example with N=5:
%
%        Zm   Ps  Zs1  P1  Zs2  P2  Zs3  P3  Zs4  P4  Zs5  P5
%    +--####--+--@@@@--|--@@@@--|--@@@@--|--@@@@--|--@@@@--+
%    |                 #        #        #        #        #    
%    0 Pin         Zp1 #    Zp2 #    Zp3 #    Zp4 #    Zp5 #
%    |                 # |->    #        #        #        #
%    +________.________|_|______|________|________|________+ 
%            x=0     x=Dx|     x=2Dx    x=3Dx    x=4Dx    x=5Dx=Lx
%			 |
%			Zeq(2)
% Ps=P(x=0); P6=P(N)=P(x=L)
%
% Related funcions TLinePs TLineVs TlineME2 
%------------------------------------------------------------------------------

% Fernando Perdigao, Coimbra, 1996 (fp@dee.uc.pt)

[Nf,Nx]=size(Zp);
[N1,N2]=size(Zs);

if (Nf~=N1), error('Size error in Zs'), end
if N2==1, Zs=Zs(:,ones(1,Nx)); end


%-- Equivalent impedances ---
%----------------------------

Zeq(:,Nx)=Zs(:,Nx)+Zp(:,Nx);	%-- at helicotrema

for k=Nx-1:-1:1,
 Zeq(:,k) = Zs(:,k) + Zp(:,k).*Zeq(:,k+1) ./ ( Zp(:,k)+Zeq(:,k+1) );
end

Zin=Zeq(:,1);		%-- coclhear input impedance (at x=0)

%--- Pressures computation with equivalent impedances
%----------------------------------------------------

Ps = Pin .* Zin./(Zin+Zm) ;				%-- Pin, Zm are scalars
P=Ps(:,ones(1,Nx)).*(cumprod((1-Zs./Zeq).').');		%-- P(k)=P(k-1)*(1-Zs(k)/Zeq(k))

%--------------------------------------------------------------