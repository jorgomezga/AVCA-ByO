  function writeHTK(fname,X,nSamples,sampPeriod,sampSize,parmKind, UNIX)
% function writeHTK(fname,X,nSamples,sampPeriod,sampSize,parmKind, UNIX)
%
% Guarda X em file de parâmetros, formato HTK (file com nome fname)
%
% sampPeriod - período de amostragem em undidades de 100ns
% paramKind  - tipo de parâmetros
% UNIX       - se UNIX~=0 : troca LSBytes com MSBytes (formato UNIX)
% X	     - matriz com (smapSize/4) linhas e nSamples colunas



%--- TESTE -------------
[i,j]=size(X);

if nSamples~=j, 
  error('writeHTK: Nº de frames diferente de nSamples')
elseif sampSize ~= i*4, 
  error('writeHTK: Nº de bytes por "amostra" diferente de sampSize')
end


if UNIX,  fp = fopen(fname,'wb','ieee-be');
else      fp = fopen(fname,'wb');
end
if fp==-1, error(['Impossivel abrir ',fname]); end

s=version; %--- na versão 5.2 do matlab, esta questão está resolvida!
	   %--- na versã0 4 não inverte bytes em formatos inteiros
if s(1)=='5' | ~UNIX,
   %--- põe header -------------
   fwrite(fp,nSamples,'int32');
   fwrite(fp,sampPeriod,'int32');
   fwrite(fp,sampSize,'int16');
   fwrite(fp,parmKind,'int16');
else
   fwrite(fp,int2byte(nSamples,4),'uchar');
   fwrite(fp,int2byte(sampPeriod,4),'uchar');
   fwrite(fp,int2byte(sampSize,2),'uchar');
   fwrite(fp,int2byte(parmKind,2),'uchar');
end

fwrite(fp,X(:),'float32');

fclose(fp);

