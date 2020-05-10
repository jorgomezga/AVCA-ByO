% argmin
% finds the arguments of the minima of matrix F(X,Y) given vector X.
% F is a vector or a matrix with the same number of rows (columns) as X
% Examples:
%   F=[3,4,2,4;  X= [1.1;         F=[2;  X=[5;
%      1,6,2,1;      2.2;            3;     6;
%      8,1,6,4]      3.3];           4]     7]
%
%    min(F,X)=[1,1,2,1]              min(F)=2;
% argmin(F,X)=[2.2,3.3,1.1,2.2]   argmin(F,X)=5
%
%function xmin=argmin(F,X)

function xmin=argmin(F,X)

[i,j]=size(F);
[k,l]=size(X);

if min(k,l)>1, error('2nd argument must be a vector'), end
if i~=k & j~=l, error('size error'), end

 
if min(i,j)==1,
 i=find(F==min(F));
 xmin=X(i(1));
elseif i==k,
 minF=min(F);
 for n=1:j,
   m=find(F(:,n)==minF(n));
   xmin(n)=X(m(1));
 end
else %--- j==l
 minF=min(F.');
 for n=1:i,
   m=find(F(n,:)==minF(n));
   xmin(n)=X(m(1));
 end
 xmin=xmin.';
end 
