  function x=miu2lin(c)
% function x=miu2lin(c)
%
% miu2lin: Converts from mu law to linear samples
%
% bits of "c" after bit invertion:
%  [s | e2 e1 e0 | m3 m2 m1 m0]
% exponent: e= e2*4+e1*2+e0
% mantissa: m= m3*8+m2*4+m1*2+m0
% signal: -1 if s=1; 1 if s=0
% x=signal*( (2^e)*(2*m+33) -33 )
%
% NOTE: MatLab function MU2LIN returns x/8192;
%       function SunAu2lin returns 4*x, 
%       where x is the vector returned by this function.
% NOTE: Maximum value in the convertion: 8031 = (2^7)*(2*15+33)-33.
%       Function MU2LIN returns always values less than 1 (max=0.98)

c=255-c; 			%inversao de bits
s=c>127;			%sinal: 0 ou 1
e = fix(c/16) - 8*s;		% expoente
m=rem(c,16);			% mantissa
x=(1-2*s).*(pow2(2*m+33,e)-33);	% x=sig*((2^e)*(2*m+33)-33)

