function [QSN,QSD,G]=tf2qs(B,A)
% Converts transfer function into quadratic sections.
%- Transforma uma funcao de transferencia dada pelos coeficientes B e A
%- em seccoes quadraticas do numerador, QSN, e do denominador, QSD.
%- G é o ganho da func. transf.
%- QSN, QSD = 	[1 a1 w12;
%		 1 a2 w22;
%		 ...
%		 0 1  wn]
%
%- correspondentes a  (s^2 + s*ak +wk2); 
%- No caso de ordem impar, a ultima linha representa (s+wn).
%- No caso de zero na origem: s*(s+wk)=(s^2 s*wk 0) -> [1 wk 0] 
%  OU:				s ->		       [0  1 0]
% [QSN,QSD,G]=tf2qs(B,A)
%versao 2
%---------------------------------------------------------------------------------


%---- normalizacao:
Gd=A(1); A=A/Gd;
G=B(1);  B=B/G;
G=G*Gd;

%-------------- zeros do numerador --------------------

z = cplxpair(roots(B));
nz=length(z);

i=1;
for k=1:2:nz-1,

 if imag(z(k))~=0,
     a1=-2*real(z(k));		% -(z1+z2)
     w02=z(k)*z(k+1);		% w02=z1*z2;
     QSN(i,:)=[1,a1,w02];
 else
   a1=-(z(k)+z(k+1));		% -(z1+z2)
   w02=z(k)*z(k+1);	% w02=z1*z2;
   QSN(i,:)=[1,a1,w02];
 end
 i=i+1;
end

if rem(nz,2),
 QSN(i,:)=[0,1,-z(nz)];
end

%----------- polos -------------

z = cplxpair(roots(A));
nz=length(z);

i=1;
for k=1:2:nz-1,

 if imag(z(k))~=0,
     a1=-2*real(z(k));		% -(z1+z2)
     w02=z(k)*z(k+1);		% w02=z1*z2;
     QSD(i,:)=[1,a1,w02];
 else
   a1=-(z(k)+z(k+1));		% -(z1+z2)
   w02=z(k)*z(k+1);	% w02=z1*z2;
   QSD(i,:)=[1,a1,w02];
 end
 i=i+1;
end

if rem(nz,2),
 QSD(i,:)=[0,1,-z(nz)];
end
