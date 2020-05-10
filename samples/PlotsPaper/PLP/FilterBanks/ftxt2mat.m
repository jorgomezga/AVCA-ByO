  function T=ftxt2mat(txtfile)
% function T=ftxt2mat(txtfile)
%
% Reads a text file into a MatLab matrix, line by line.
% (usefull to insert/delete text by columns (! tabs))
% 
% function T=ftxt2mat(txtfile);
%
% To convert T again to a text file, use: mat2ftxt(txtname,T)
%

f=fopen(txtfile);
if f==-1, error(['Cannot open ',txtfile]); end

s=fread(f,'char'); fclose(f);

i=[find(s==13),find(s==10)]';	%--- <CR> and/or <LF>
D=size(i,1);
if D>1, i=max(i); end
N=length(i);

n=length(s);
if n==i(length(i)), i=[0,i];
else i=[0,i,n+D]; N=N+1;
end
len=diff(i)-D;	%-- excludes <CR> and/or >LF>
MaxL=max(len);

T = ones(N,MaxL)*' ';

for k=1:N,
 if len(k),
  T(k,1:len(k))  =  s(i(k)+1:i(k+1)-D)';
 end
end

T=setstr(T);

