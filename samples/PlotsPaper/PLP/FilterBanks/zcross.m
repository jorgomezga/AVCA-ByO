% ZCROSS 
% Given the vectors f and x, ZCROSS finds the zero crossings of f(x) and
% returns the linearly interpolated values xz from the values of x 
% for which f crosses zero or is zero.
%
% function xroots=zcross(f,x)
%
%
% Examples:
%   f=[4 9 16 25 36]
%   x=[2 3  4  5  6]   ,  f(x)=x^2
%
%   xz = zcross(f-5,x) returns xz=2.2
%
% t=0:0.5:7;
% t0=zcross(sin(t),t) returns t0=[0;3.1434;6.2825]
%

function xz=zcross(f,x)

N=length(f);
if length(x)~=N, error('f and x have different sizes'); end

f=f(:); x=x(:);
fa=f(1:N-1);
fp=f(2:N);


iz=find(f==0);
i=[find(fa<0 & fp>0) , find(fa>0 & fp<0)];
if length(i),
   xz=x(i)-f(i).*(x(i+1)-x(i))./(f(i+1)-f(i));
end
xz=sort([x(iz);xz]);
