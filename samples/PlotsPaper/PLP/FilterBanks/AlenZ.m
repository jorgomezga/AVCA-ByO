function [Zp,Ht]=AllenZ(f,x);
% function [Zp,Ht]=AllenZ(f,x);
% Computes cochlear partition impedance for Allen's model.
% Returns also Ht=VC/VBM the relationship between cilar and 
% basilar membrane velocities.
% 
% Ref.:J. Allen, "Cochlear Micromechanics-A Physical Model of Transduction"
% JASA 68(6),pp.1660-1669,1980. Also JASA 89(1),pp.287-309, 1991.
%


%--- model constants ----
M=0.04;
Mb=M/1.2;
Mt=0.2*Mb;
L=2.5;
%--------------------------


ux=ones(size(x));
uf=ones(size(f));
s=j*2*pi*f;
invs = 1.0./s;

G=0.5*exp(1-x/L);
wCF=2*pi*456*(10 .^(2.1*(1-x/L))-0.9999);
wp=1.3*wCF;
wz=0.65*wCF;
dp=2*0.3;
dz=2*0.5;

wp2=wp.^2;
wz2=wz.^2;


%---- K/M = wCF^2 => K = Kb + (G^2)*Kc*(wz/wp)^2
%
Kc = Mt*(wp2-wz2);
Rc = Mt*(dp*wp-dz*wz);
Kt = Mt*wz2;
Rt = Mt*dz*wz;
Kb = M*(wCF.^2) - (G.^2).*Kc.*(wz2./wp2);
Kb=3*Kb; %--- as in Puria, 91.


Zt = Mt*(s*ux) + uf*Rt + invs*Kt;
Zc = uf*Rc + invs*Kc;
Ht=(uf*G).*(Zt./(Zc+Zt));

ZBM = M*(s*ux) + invs*Kb;
Zp  = ZBM + (uf*G).*Zc.*Ht;

