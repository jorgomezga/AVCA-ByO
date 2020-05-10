% function mod = db(H,Hmin)
%
% dB(H)=20*log10(abs(H))
% Hmin: log floor (default=1e-10). Ex: db(0,1e-5)=-100

function mod = db(H,Hmin)

if nargin<2, Hmin=1e-10; end
mod = abs(H);
i=find(mod<Hmin);
l=length(i);
if l==1, mod(i) = Hmin;
elseif (l>1), mod(i) = Hmin * ones(1,l); end
mod = 20*log10(mod);
