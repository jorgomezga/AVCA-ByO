  function [X,nSamples,sampPeriod,sampSize,parmKind]=readHTK(fname,UNIX)
% function [X,nSamples,sampPeriod,sampSize,parmKind]=readHTK(fname,UNIX);
%
% Lê file "fname" para matriz X  - file de parâmetros, formato HTK
%
% Toma header de forma a obter os seguintes parâmetros:
% nSamples   - numero de frames
% sampPeriod - período de amostragem em undidades de 100ns
% sampSize   - Dimensão da frame * 4 (num. de bytes de 1 frame)
% paramKind  - tipo de parâmetros
%
% UNIX       - 0: formato PC; 1: formato UNIX (ordem dos bytes trocada)
% X	     - matriz com (sampSize/4) linhas e nSamples colunas


if UNIX, fp = fopen(fname,'rb','ieee-be');
else      fp = fopen(fname,'rb');
end
if fp==-1, error(['Impossivel abrir ',fname]); end

s=version; %--- na versão 5.2 do matlab, esta questão está resolvida!

if s(1)=='5' | ~UNIX,
   %--- traz header -------------
   nSamples  = fread(fp,1,'int32');
   sampPeriod= fread(fp,1,'int32');
   sampSize  = fread(fp,1,'int16');
   parmKind = fread(fp,1,'int16');
else
   %---- na versão 4 o big-endian só funciona com floats
   h=fread(fp,12,'uchar');
   nSamples   = hex2dec(sprintf('%x%x%x%x',h(1:4)));
   sampPeriod = hex2dec(sprintf('%x%x%x%x',h(5:8)));
   sampSize   = hex2dec(sprintf('%x%x',h(9:10)));
   parmKind  = hex2dec(sprintf('%x%x',h(11:12)));
end

X=fread(fp,[sampSize/4,nSamples],'float32');
%---- verifica tamanho -------
[i,j]=size(X);
if i~=sampSize/4 | j~=nSamples,
 fclose(fp); error('ficheiro de parâmetros incorrecto');
end

fclose(fp);

%teste
%return
%C=ones(3,1);
%C=[C*1.1111,C*2.2222,C*3.3333,C*4.4444,C*5.5555];
%nSamples=5;
%sampPeriod=100000;
%sampSize=3*4;
%parmKind=6;
%writeHTK('lixo.par',C,nSamples,sampPeriod,sampSize,parmKind, 1)
%[X,nSamples,sampPeriod,sampSize,parmKind]=readHTK('lixo.par',1)

