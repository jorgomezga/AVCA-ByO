 function [P,Ps,Zin,Zeq]=TLineME2(Pin,Zm,Zs,Zp)
%function [P,Ps,Zin,Zeq]=TLineME2(Pin,Zm,Zs,Zp)
%
% Transmission Line equivalent model of the cochlea with middle ear connected.
% Computes cochlear pressures/forces and the equivalente impedances.
%
% Input: series and parallel distributed impedances:
%	 series impedance: Zs: Zs(Nf,1) or Zs(Nf,Nx), Nf: no. of points in f.
%	 parallel impedances: Zp(Nf,Nx),              Nx: no. of points in x.
%	 Zm:	Middle ear impedance referred to stapes velocity
%	 Pin:	Input pressure/force "	"	"
% Output: P:	Pressures/Forces  P(1)..P(Nx)		(Nf x Nx)
%	  Ps:	Fluid pressure/force at stapes		(Nf x 1)
%	  Zin:	Input impedance of the cochlea: 	(Nf x 1)
%	  Zeq:	Equivalent impedances			(Nf x Nx)
%
% Note: See also TLineME.m (differs in the definition of Zeq)
%------------------------------------------------------------------------------
%
% Nx=5 sections:
%
%        Zm   Ps  Zs1  P1  Zs2  P2  Zs3  P3  Zs4  P4  Zs5  P5
%     --####--+--@@@@--|--@@@@--|--@@@@--|--@@@@--|--@@@@--+
%    +                 #        #        #        #        #    
%   Pin            Zp1 #    Zp2 #    Zp3 #    Zp4 #    Zp5 #
%    -             |-> #        #        #        #        #
%     ________.____|___|________|________|________|________+ 
%            x=0   |  x=Dx      x=2Dx    x=3Dx    x=4Dx    x=5Dx=Lx
%		Zeq(1)
%
% Notes:
% Ps=P(x=0); P(Nx)=P(x=L) (at helicotrema)
% Zin=Zs1+Zeq(1) - input impedance of the cochlea (specific if pressures)
% With fluid velocities, Vx(k)=(P(k-1)-P(k))/Zs(k), then Zeq(k)=P(k)/Vx(k)
%
% Related funcions TLinePs TLineVs TLineME TlineME2 
%------------------------------------------------------------------------------

% Fernando Perdigao, Coimbra, 1996 (fp@dee.uc.pt)

function [P,Ps,Zin,Zeq]=TLineM(Pin,Zm,Zs,Zp)


[Nf,Nx]=size(Zp);
[N1,N2]=size(Zs);

if (Nf~=N1), error('Size error in Zs'), end
if N2==1, Zs=Zs(:,ones(1,Nx)); end

[N1,N2]=size(Pin);
if length(Pin)>1,
 if ((N1~=Nf)|(N2>1)), error('Size error in Pin'), end;
 Pin=Pin(:,ones(1,Nx));
end

%-- Equivalent impedances ---
%----------------------------

Zeq(:,Nx) = Zp(:,Nx);	%-- at helicotrema

for k=Nx-1:-1:1,
 Zeq(:,k) = Zp(:,k).*( Zs(:,k+1)+Zeq(:,k+1) )./( Zp(:,k)+Zs(:,k+1)+Zeq(:,k+1) );
end

Zin=Zs(:,1)+Zeq(:,1);		%-- coclhear input impedance (at x=0)

%--- Pressures computation with equivalent impedances
%----------------------------------------------------

Ps = Pin .* Zin./(Zin+Zm) ;			%-- Pin and Zm have the same dimension

P=Ps./(cumprod((1+Zs./Zeq).').');	%-- P(k)=P(k-1)*Zeq(k)/(Zeq(k)+Zs(k))

