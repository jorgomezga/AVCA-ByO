  function [A,B,h0,Vcenter]=MedConst(Th,fspo,fsat,NewConst)
% function [A,B,h0,Vcenter]=MedConst(Th,fspo,fsat,NewConst);
%
% Computs the Meddis' IHC model constants A, B and h0,
% in order to have a model with absolut threshold Th (dB); 
% spontaneous firing rate fspo and saturation firing rate fsat.
%
% The other model constants are left unchanged (default values) 
% unless new constants are given in vector NewConst.
% NewConst=[g0,y0,x0,l0,r0]
% default constants: g0=2000; y0=5.05; x0=66.31; l0=2500; r0=6580.
%
% Returns also Vcenter, the sinusoidal amplitude at 1kHz that 
% produces a firing rate f=(fsat-fspo)/2.
% Displays model characteristics and adaptation time constants.
% Ref: "Implementation Details of a Computational Model of the
%       Inner Hair Cell/Auditory-nerve Synapse", 
%       JASA, 87(4) Apr. 1990, pp.1813-1816.

% F. Perdigao, Coimbra, 1997
% fp@it.uc.pt


A = 10^(Th/20);

if (fsat<fspo) error('fsat must be higher than fspo'); end

if nargin<4,	%--- Default Constants
 g0 = 2000;
 y0 = 5.05;
 x0 = 66.31;
 l0 = 2500;
 r0 = 6580;
else		%--- New Constants
 g0=NewConst(1);
 y0=NewConst(2);
 x0=NewConst(3);
 l0=NewConst(4);
 r0=NewConst(5);
end


D= y0/l0;
E= (l0+r0)*y0/l0;
%---------------------------------------------------


cmax= D*g0/2/(g0/2+E);
h0=fsat/cmax;

%--- fspo = h0*D*(g0*F) / [g0*F + E], F=A/(A+B)
 
F=fspo*E/g0/(D*h0-fspo);
B=A/F-A;


if any([A,B,h0]<0), error('constants less than zero'), end


%--- RATE-INTENSITY PLOT --------

Vdb= linspace(Th-20,Th+80,500);
V=10.0.^(Vdb/20);
kV = mean_k(V,A,B,g0);
fV = h0*D*kV./(E+kV);

plot(Vdb,fV)
title('Rate-Intensity curve'); xlabel('V [dB]'); ylabel('rate');



%------- Dynamic range in 5% to 95% of firing rate

fcenter = (fsat+fspo)/2;
i=find_ind(fV,fcenter);
Vcenter=Vdb(i);

v1=zcross(fV-1.05*fspo,Vdb);
v2=zcross(fV-0.95*fsat,Vdb);

disp('----------------------------------------------------------');
disp('Meddis'' IHC/Synapse model');
disp(['A=',num2str(A,3),' B=',num2str(B,3),' h0=', num2str(h0,3)]);
disp(['Dynamic range (5% to 95% of absolute values): ', num2str(v2-v1,3), ' dB']);
disp(['Rate threshold: ', num2str(v1,3), ' dB']);
disp(['Sat. threshold: ', num2str(v2,3), ' dB']);
disp(' ');

%--- Time constants ---------------
A1=[1,l0+r0];
B2=[1,(l0+r0+x0),l0*x0];
A2=conv([1,y0],[1,x0]);
C=A+B;
kmin=g0*A/C; 
kmax=g0/2;

A3=conv(A2,A1)+kmin*[0,B2];
B3=kmin*A2;
P=roots(A3);
tau_min=-1 ./P;

A3=conv(A2,A1)+kmax*[0,B2];
B3=kmax*A2;
P=roots(A3);
tau_max=-1 ./P;

disp('Time Constants for a burst with instant. rise time:');
fprintf('Short-Term adaptation: tau1=%5.1fms (V=0) ... %5.1fms (V=Inf)\n',...
        tau_min(2)*1000,tau_max(2)*1000);
fprintf('Rapid adaptation     : tau1=%5.1fms (V=0) ... %5.1fms (V=Inf)\n\n',...
        tau_min(3)*1000,tau_max(3)*1000);

%--- aprox. values directly from constants: ------------
t1min=1000/(g0/2+x0+y0);
t1max=1000/(g0*A/C+x0+y0);
t2min=(l0+r0)/x0/l0*1000;
t2max=t2min/t1max/(g0*A/C+E)*1000-t1max;

disp('Aprox. values:');
fprintf('tau1=%5.1fms ... %5.1fms\n',t1max,t1min); 
fprintf('tau2=%5.1fms ... %5.1fms\n\n',t2max,t2min); 
disp('----------------------------------------------------------');

