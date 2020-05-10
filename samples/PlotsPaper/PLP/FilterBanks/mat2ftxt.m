  function mat2ftxt(txtfile,T);
% function mat2ftxt(txtfile,T);
%
% Writes a text file from a MatLab matrix.
% 
% see ftxt2mat
%

f=fopen(txtfile,'wt');
if f==-1, error(['Cannot open ',txtfile]); end

for k=1:size(T,1),
  fprintf(f,'%s\n',deblank(T(k,:)));
end

fclose(f);

