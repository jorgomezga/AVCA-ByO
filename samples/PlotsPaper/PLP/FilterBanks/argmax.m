% argmax
% finds the arguments of the maxima of matrix F(X,Y) given vector X.
% F is a vector or a matrix with the same number of rows (columns) as X
% Examples:
%   F=[3,4,2,4;  X= [1.1;         F=[2;  X=[5;
%      1,6,2,1;      2.2;            3;     6;
%      8,1,6,4]      3.3];           4]     7]
%
%      max(F)=[8,6,6,4]               max(F)=4;
% argmax(F,X)=[3.3,2.2,3.3,1.1]   argmax(F,X)=7
%
%function xmax=argmax(F,X)

function xmax=argmax(F,X)

[i,j]=size(F);
[k,l]=size(X);


if min(k,l)>1, error('2nd argument must be a vector'), end
if i~=k & j~=l, error('size error'), end

 
if min(i,j)==1,
 i=find(F==max(F));
 xmax=X(i(1));
elseif i==k,
 maxF=max(F);
 for n=1:j,
   m=find(F(:,n)==maxF(n));
   xmax(n)=X(m(1));
 end
else %--- j==l
 maxF=max(F.');
 for n=1:i,
   m=find(F(n,:)==maxF(n));
   xmax(n)=X(m(1));
 end
 xmax=xmax.';
end 
