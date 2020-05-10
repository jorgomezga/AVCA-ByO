% J. ROSOWSKI EXTERNAL AND MIDDLE-EAR MODEL
% Ref: "Models of External and Middle-Ear Function"
%      in "Auditory Computation", SHAR-6, Springer, 1993
%
% function [Zeq,Peq,PC,PS,PT,PEC,PEX] = rosowski(f,PPW)
%
% PPW - Plane wave pressure
% f   - frequency values
% Outputs:
% PEX - external pressure (at pinna)
% PEC - pressure at ear canal
% PT  - pressure at timpanic membrane
% PC  - pressure at the base of the cochlea
%
% Peq - Thevenin pressure,
% Zeq - Thevenin impedance, according to the following model:
%
%       Zeq   --> US         US: stapes volume velocity 
%   +--#####----+------+     PC: pressure at the base of the cochlea
%   |                  |     ZC: cochlear input impedance
%   |                  #         
%  (+)Peq       PC     # ZC
%   |                  # 
%   |                  |
%   ------------+------+
%
% NOTE: uses Bauer's 67 model for sound diffraction and 
%       scattering by the head.
%       (see Giguere, & Woodland, JASA (95)-1,1994, pp.332)
% NOTE: several values of original model was changed, namely
%       the cochlear input impedance.

% F. Perdigao, IT-Coimbra, Dec.1996


function [Zeq,Peq,PC,PS,PT,PEC,PEX,ZiT] = rosowski(f,PPW)


if nargin <2, PPW=1; end

%% PPW=1; f=logspace(1,5,500);

s=j*2*pi*f;
w=2*pi*f;
c=343;		%--- propagation velocity of sound in air
r0=1.21;	%--- air density
A0=4.4e-4;	%--- concha horn area (wide end)
Al=A0/10;	%--- concha horn area (narrow end)
lH=0.01;	%--- horn length
lT=0.02;	%--- tube length
as=0.1;		%--- radius of head sphere

%=============================================
%--- 1st part: The Equivalent Pressure Source,
%    Diffraction and Scattering by the head
%=============================================

%--- Simplified model of Bauer, 67 
%--- (see Giguere, & Woodland, JASA (95)-1,1994, pp.332)
%-------------------------------------------------------

Zh=s*r0/(2*pi*as);
Rh=r0*c/(pi*(as^2));
Zr=s*0.7*r0/sqrt(pi*A0);
Rr=r0*c/A0;

GS=(Rh+2*Zh)./(Rh+Zh);		%--- 6dB gain at high freqs.
ZR= zparal(Rh,Zh)+zparal(Rr,Zr);%--- Thevenin source impedance

%=============================================
%--- 2nd part: Concha and External Canal,
%=============================================

%--- HORN

alpha=0;
k=w/c-j*alpha;
a=0.5*log(Al/A0)/lH;
b=sqrt(k.^2 - a^2);
k1=-a-j*b;
k2=-a+j*b;
AH=exp(a*lH)*(k2.*exp(j*b*lH)-k1.*exp(-j*b*lH))./(2*j*b);
BH=exp(a*lH)*s*r0.*sin(b*lH)./(Al*b);
CH=-exp(-a*lH)*Al*(a^2+b.^2).*sin(b*lH)./(s*r0.*b);
DH=exp(-a*lH)*(k2.*exp(-j*b*lH)-k1.*exp(j*b*lH))./(2*j*b);

%--- TUBE

k=w/c-j*alpha;	
z0=r0*c/Al;
AT=cos(k*lT);
BT=j*z0*sin(k*lT);
CT=j/z0*sin(k*lT);
DT=AT;

A=AH.*AT+BH.*CT;
B=AH.*BT+BH.*DT;
C=CH.*AT+DH.*CT;
D=CH.*BT+DH.*DT;	%--- A.*D-B.*C = 1 (reciprocity)

%======================================================
%--- 3rd part: Thevenin Equivalent at timpanic membrane
%======================================================
%
% GS*PPW = ZR*UEX + PEX;
%        = ZR*(C*PT+D*UT) + (A*PT+B*UT)
%        = PT*(ZR*C+A) + UT*(ZR*D+B)

Pem = PPW.*GS./(ZR.*C+A);
Zem = (ZR.*D+B)./(ZR.*C+A);
semilogx(f,db(Pem))

%======================================================
%--- 4th part: Thevenin Equivalent at stapes
%======================================================

%--- Middle-Ear impedances

ATM=60e-6;
AFP=3.2e-6;
rMI=1.3;
CMC=3.9e-11;
CTC=4e-12;
RA=6e6;
LA=100; LT1=750; LT=6.6e3;
CT=3e-12; CT2=1.3e-11;
RT2=1.2e7;
%RT=2e6; 	%--- original
RT=3e7;		%--- altered


CTS=1.1e-3;	%--- original
%CTS=5e-3;	%--- altered
RTS=4.3e-2;
RMI=7.2e-3; 
LMI=7.9e-6;
%CMI=Inf;
LS=3e-6;
RJ=3.6; 
CJ=4.9e-4;
CAL=9.4e-15;	%--- original
CAL=5e-15;	%--- altered

%---- COCHLEAR INPUT IMPEDANCE
%=============================
%---- cochlear input impedance (wrong!)
%--RC=2e10; LC=2.4e6;
%--ZC= RC+s*LC;

%-- from Lynch model
%R0C=0.28e11;RC=1.2e11;LC=2.25e8;CRW=1e-13;MV=22e5;
%CAL=0.37e-14;MS=3.3e5;RAL=0.2e11; 
%ZC=zparal(RC,R0C+s*LC)+1.0./(s*CRW)+s*MV;ZAL=RAL+s*MS+1.0./(s*CAL);

%-- from Neely's model
%R0C=2.2e10;RC=1.55e11;LC=8e7;

R0C=5e10;RC=1e11;LC=8e7; 		%-- used impedance!

CRW=1e-13;
ZC=zparal(RC,R0C+s*LC) +1.0./(s*CRW);
%ZAL=1.0./(s*CAL);
ZAL=1.0./(s*CAL)+1e10;		%------- with added resistor

ZTC=zparal(1.0./(s*CTC) , s*LA+RA+1.0./(s*CMC));
ZT = s*LT1 + zparal(1.0./(s*CT)+RT+s*LT , 1.0./(s*CT2)+RT2 );
ZTS=1.0./(s*CTS)+RTS;
%ZMI=RMI+s*LMI+1.0./(s*CMI);
ZMI=RMI+s*LMI;			%--- with CMI=Inf
ZJ=RJ+1.0./(s*CJ);


%---- Reduction to transformer secundaries  --------
%===================================================

%---- 1st tranformer ---------
Peq=ATM*Pem;
Zeq=(ATM^2)*(Zem + ZT + ZTC);

%---- 2nd transformer --------
Peq = rMI*Peq.*ZTS./(ZTS+Zeq);
Zeq = (rMI^2)*(ZMI + zparal(Zeq,ZTS));

%---- 3rd transformer --------
Peq = (1/AFP)*Peq.*ZJ./(ZJ+Zeq);
Zeq = (1/(AFP^2))*(s*LS + zparal(Zeq,ZJ));

%---- Series Equivalent Impedance seen at cochlea ---
Zeq = Zeq + ZAL + 1.0./(s*CRW);


%--- PC/PPW and PS/PPW --------
%==============================

PC=Peq.*ZC./(ZC+Zeq+ZAL);
PS=Peq.*(ZC+ZAL)./(ZC+Zeq+ZAL);


%---- ZiT - Impedance at Timpanus
%=================================

ZiT = (ZAL+ZC)*(AFP^2); 
ZiT = zparal(s*LS + ZiT , ZJ);
ZiT = ZiT/(rMI^2);
ZiT = zparal(ZMI+ZiT,ZTS);
ZiT = ZiT/(ATM^2);
ZiT = ZT+ZiT+ZTC;

%----- PT/PPW -----------
%========================

PT = PPW.*GS./(ZR.*(C+D./ZiT)+A+B./ZiT);

%----- PEC and PEX ----------
%============================

PEC = PT.*(AT+BT./ZiT);
PEX = PT.*(A+B./ZiT);


%ZiEX = (A+B./ZiT)./(C+D./ZiT);
%ZiEC = (AT+BT./ZiT)./(CT+DT./ZiT);
%ZoEX=(B+D.*ZR)./(A+C.*ZR);

semilogx(f,db(PC),f,db(PS),f,db(PT),f,db(PEC),f,db(PEX))

return



%------- Model for Zeq: ------------

Zeqmod=s*1.4e6+ 2e10 + 2.3e14./s;
semilogx(f,db(Zeq+ZAL),f,db(Zeqmod))

%-- more detailed:
b2 =[2.8939e+005,4.6250e+011,3.4374e+016,3.2750e+021,1.4567e+026,4.2793e+030,...
     8.7703e+034,1.2144e+039,9.6269e+042,4.6572e+046,5.1613e+049];
a2 =[1.0000e+000,3.3878e+005,2.3356e+010,1.9967e+015,7.2365e+019,1.6808e+024,...
     2.6765e+028,1.7384e+032,2.1299e+035,0];
h=freqs(b2,a2,w);
semilogx(f,db(Zeq),f,db(h))

%-------------------------------------

%-----Model for Peq -----------------------

[b,a]=invfreqs(Peq,w,6,6,[ones(1,400),zeros(1,100)],10);
b =[1.0683e+000 -7.1108e+005  4.8589e+010 -1.9088e+015 -4.6478e+019  3.2648e+024 7.8505e+027];
a =[1.0000e+000  6.5762e+004  6.7790e+009  2.8729e+014  6.2436e+018  1.0478e+023 5.0628e+026];

h=freqs(b,a,w);
semilogx(f,db(Peq),f,db(h))
%-------------------------------------------

