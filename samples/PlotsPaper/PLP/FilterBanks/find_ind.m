%FIND_IND
%  Given the vector f and constant (or vector) c,
%  FIND_IND finds the INDEX of the first occurrences 
%   of f(i)~=c, in the sense:
%       |f(i)-c| = min(|f-c|) 
%
% Example:
%   f=[10 20 30 40 50 60]
%   c=[22 29]
% returns: index=[2 3]
% and f(index)=[20 30].
%
% function index = find_ind(f,c)
%


function index = find_ind(f,c)

nc=length(c);

for k=1:nc,
 i=find(abs(f-c(k)) == min(abs(f-c(k))));
 index(k)=i(1);
end