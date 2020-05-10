% function [Zin,Zeq]=coc_zin(Zs,Zp);
%
% Computs the input impedance of a transmission line cochlear model:
%
%            Zs(1)   Zs(2)    Zs(3)       Zs(N) 
%         --@@@@--|--@@@@--|--@@@@--|-...--@@@@--|
%                 #        #        #            #   
%   Zin -->  Zp(1)#   Zp(2)#   Zp(3)# ...        # Zp(N)
%                 #        #        #            #
%         ________|________|________| ...________| 
%
% Input:  Partition impedances Zp(Nf,Nx), Zs: Zs(Nf,1) or Zs(Nf,Nx)
% Output: Zin=Zeq(:,1) (Nf x 1) and  Zeq (Nf x Nx)
%-------------------------------------------------------------------


function [Zin,Zeq]=coc_zin(Zs,Zp);

%-- Equivalent impedances ---
%----------------------------

[Nf,Nx]=size(Zp);
[N1,N2]=size(Zs);

if (Nf~=N1), error('Size error in Zs or Zp'), end
if N2==1, Zs=Zs(:,ones(1,Nx)); end


Zeq(:,Nx)=Zs(:,Nx)+Zp(:,Nx);
for k=Nx-1:-1:1,
 Zeq(:,k) = Zs(:,k) + Zp(:,k).*Zeq(:,k+1) ./ ( Zp(:,k)+Zeq(:,k+1) );
end

Zin=Zeq(:,1);
