  function x=readpcm(filename)
% function x=readpcm(filename)
%
% Reads file with raw data PCM mu law signal into x (linear samples).
 
  fp=fopen(filename,'r');
  if fp==-1, error(['Cannot open ',filename]); end
  x=fread(fp,inf,'uchar');
  fclose(fp);
  x=miu2lin(x);


