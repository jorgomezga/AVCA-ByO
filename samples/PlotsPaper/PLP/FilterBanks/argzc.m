% Finds the zero crossings of f(x) and
% returns the x arguments where f crosses zero,
% from negative to positive and from positive to negative values
%
% function [xa,xd]=argzc(f,x)
% xa: x for ascendent crossings
% xd: x for descendent crossings
%
% if x is the indexes of f, use: [i,j]=argzc(f)


function [xa,xd]=argzc(f,x)

N=length(f);
f1=f(1:N-1);
f2=f(2:N);

i=find((f1<=0)&(f2>0));
j=find((f1>=0)&(f2<0));

if nargin >1
 xa=x(i);
 xd=x(j);
else
 xa=i;
 xd=j;
end