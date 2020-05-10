% lin2miu - converts samples to mu law.
%
% |x|=(2^e)*(2*m+33)-33  with n (8bits): [s e2 e1 e0 m3 m2 m1 m0]
% n= s*128 + e*16 + m
% u=255-n (bit-inverted sample)

function n=lin2miu(x)

s=x<0;
x=min(abs(x),8031);
a=x+33;			%-- a=(2^e)*(2*m+33)
b=fix(a/33);		%-- b=(2^e)*(1+y)  com y<1
e=fix(log2(b));
m=fix(a./pow2(e+1))-16;  %-- mantissa. a/(2^(e+1))=(m+33/2)
n=255-(128*s+16*e+m);

