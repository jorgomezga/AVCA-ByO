  function [A,B,GHW,uf,fmax]=SenConst(VC,     fspo,fsat,Vdb)
% function [A,B,GHW,uf,fmax]=SenConst(Vcenter,fspo,fsat,Vdb);
%
% Evaluates the Seneff's model constants, A, B and GHW
% in order to obtain a fiber model with spontaneous rate (fspo)
% and saturated rate (fsat) and a rate-intensity curve with center
%         fcenter=(fsat+fspo)/2; V=Vcenter.
% The rate-threshold is aprox.: Vcenter -25dB. (Vcenter in dB).
%
% The other model constants are not changed in order to have aprox. 
% the same adaptation time constants.
%  
% Vdb is an optional vector with amplitude values for the plot
% of the rate-intensity curve. If "uf" and/or "fmax" are specified,
% then the model is evaluated for the Vdb values and returns
% uf, the rate curve, and fmax, the pick of the mean rate (on 1ms basis)
% for 1kHz input, 2.5ms rise time, and with the model initially at 
% rest conditions.
  
% Fernando Perdigão, DEEUC/IT, Coimbra 1997
 

if (fsat<fspo) error('fsat must be higher than fspo'); end
if (fsat >400) error('fsat must be less than 400 for the same KAGC=0.002'); end
if (fspo<=0) error ('fspo <=0'); end

wb=1000/120;
wa=1000/15 - wb;
KAGC=0.002;
c=wa*wb/(wa+wb);

zspo = fspo/(1-KAGC*fspo);	%--- F=Z/(1+KAGC*Z); Z=F/(1-KAGC*F);
GHW  = zspo/c;			%--- zspo=GHW*[wa*wb/(wa+wb)]
zsat = fsat/(1-1.01*KAGC*fsat);
A    = (zsat/GHW-6.28)/9.87;	%--- zsat=GHW*(9.87*A+6.28) ~= 
				%--- ~= GHW*(1+A*pi/2)*wa*wb/(wa+2*wb);

%--------------------------------------------------------------------
%--- Constant B from the center of the rate-intensity curve ---------
%--------------------------------------------------------------------

%-- empirical constants for wa, wb and KAGC fixed, B=1:----------------------
%-- 2 parameters to caracterize mean(f) as zspo+(zsat-zspo)*V^n/(V^n+Vc^n)
%-- as a function of A.
n =1.146+0.0053*exp(-0.0181*A)+0.3485*exp(-0.8280*A)+0.0252*exp(-0.1267*A);
Vc=1.390+0.0260*exp(-0.0243*A)+0.4234*exp(-0.5287*A)+0.1028*exp(-0.1285*A);
%----------------------------------------------------------------------------

fc=(fsat+fspo)/2;		%--- rate at center
zfc=fc/(1-fc*KAGC*1.01);	%--- equivalent mean(z)
c=(zfc-zspo)/(zsat-zspo);
VC1=db( c*(Vc^n)/(1-c) )/n;	%--- Amplitude V (in dB) at the center for B=1
dbB=VC1-VC;			%--- B for the desired center
B=10^(dbB/20);

if any([A,B,GHW]<0), error('constants less than zero'), end
if nargin>3,
 NV=length(Vdb);
else 
 NV=200;
 Vdb=linspace(VC-60,VC+50,NV);
end

V=10.0.^(Vdb/20);

%---------- aproximate curves for mean(z) and mean(f) -----------------
v=V*B;
uzmod = zspo+(zsat-zspo)*(v.^n)./(v.^n+Vc^n);
ufmod = uzmod./(1+uzmod*KAGC*1.01);
%----------------------------------------------------------------------

if nargout > 3,	%--- compute the empirical rate curve
		%#########################################
		% NOTE: Requires a lot of time to evaluate
		%#########################################
fs = 16000;
T = 1/fs;
t = (0:T:0.4-T);		%--- 400 ms burst
N=length(t);
env = 1 - exp(-t/0.0025);	%--- burst envelope w/ 2.5ms rise time
f0=1000;			%--- input frequency
N0=fs/f0;			%--- no. samples per period
N1ms=N/N0; 			%--- no. of 1ms periods

burst = sin(2*pi*f0*t);		%--- sinusoidal burst at f0

for i=1:NV,
i
 x=V(i)*burst.*env;
 [f,z]=Seneff(x,fs, A, B, GHW);
 fav = mean(reshape(f,N0,N1ms));	%--- mean(c) em períodos inteiros de 1ms
 zav = mean(reshape(z,N0,N1ms));
 uf(i)=fav(N1ms);
 uz(i)=zav(N1ms);
 fmax(i)=max(fav);
end %--- i

plot(VdB, uf)

else
plot(Vdb,ufmod,'m')
end		%if nargout>3


v1=zcross(ufmod-1.05*fspo,Vdb);
v2=zcross(ufmod-0.95*fsat,Vdb);

disp('Seneff''s IHC/Synapse model:');
disp(['A=',num2str(A,3),' B=',num2str(B,3),' GHW=', num2str(GHW,3)]);
disp(['Dynamic range (5% to 95% of absolute values): ', num2str(v2-v1,3), ' dB']);
disp(['Rate threshold: ', num2str(v1,3), ' dB']);
disp(['Sat. threshold: ', num2str(v2,3), ' dB']);



