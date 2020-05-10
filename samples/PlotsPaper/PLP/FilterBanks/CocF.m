% CocF
% Cochlear (Functional) Model
% Computes an auditory model for 8KHz sampling rate
% Uses 35 sections of cascade filters as well as parallel filters
% Output: Filter coeficients:
% Bm, Am, Bt,At for cascade filters
% Bs As for parallel filters
%
% Output: Filter coeficients:
%         Bm, Am, Bt,At for cascade filters
%	  Bs As for tap filters
%
%         +------+     +-------+                  +--------+
%  x ---->|  Hm  |---->| Ht(1) |--+--> ... ------>| Ht(35) |---+
%         +------+     +-------+  |               +--------+   |
%                             +-------+                    +--------+
%                             | Hs(1) |       ...          | Hs(35) |
%                             +-------+                    +--------+
%                                 |                            |
%                                y(1)         ...            y(35)
%
% Ref:
% F. Perdigao, L.V. de Sa', 
%"A Cochlear Model for Speech Processing", RecPad’95, Aveiro, 1995



fs=8000;                     % sampling frequency
fN=fs/2;
NS=35;                       % no. of sections
NPOINTS=512;                 % points for freq. responses

%------- defining f and W ----------------------
f  = logspace(2,log10(fN),NPOINTS);
W  = 2*pi*f/fs ;
fapex = 300;		%-- CF of last section
fbase = 3400;		%-- CF of first section


%--- defining the interest range of cochlear positions (0<x<1) --------------
xapex =  log10(fapex/160 + 0.8)/2.1;	% cochlear place equivalence (Liberman)
xbase =  log10(fbase/160 + 0.8)/2.1;	% "base" or "stapes"

%-------- cochlear place, x, is linearly spaced -----

x = linspace(xbase,xapex,NS);


%-------- Liberman's freq. to place function (see Greenwood,90) ------

fCF  = 160*(10 .^ (2.1*x) -0.8);
WCF  = 2*pi*fCF/fs;		%--- center freq.s (when |Yk(W)| is maximum)
fiCF = tan(WCF/2);		%--- bilinear warping freq.: fiCF=wCF/(2*fs)=tan(WCF/2)


%---- normalized x (to define parameter variations) --------------

%%fnorm = (fCF-fapex)/(fbase-fapex);
xnorm = (x-xapex)/(xbase-xapex);


%---- fCTAP: freq. of cascade tap responses's maxima. ----
%---- Pk(W) is the freq. response of the cascade at tap k
%---- Yk(W) is the freq. response of the parallel output k of the model
%---- The max of |Pk(W)| is bellow the max of |Yk(W)| in order to obtain
%---- sharper responses above CF.
%---------------------------------------------------------------------

fCTAP=fCF-15;		%---- when |Pk(W)| is maximum
WCTAP=2*pi*fCTAP/fs;
fiCTAP=tan(WCTAP/2);

%--- Filters' Quality factors -------------
%--- The pole Q's will be reestimated; the zero Q's remain as defined here

Qps = 9 + (5-9)*xnorm;		%--- from 5 to 9
Qps(1)=4; 			%--- change the first 2 Q's
Qps(2)=4;

%Qzs = 4 + (4-4)*xnorm;
Qzs = 1 + (4-1)*xnorm;

Qpt=1.8+(3-1.8)*xnorm;		%--- from 3 to 1.8
Qpt(1)=1.6;

%%Qzt=5+(8-5)*xnorm;
Qzt=2+(6-2)*xnorm;

%--- equivalent damping (or two times damping) factors ----------------- 

dpt = 1 ./Qpt;  dzt = 1 ./Qzt;
dps = 1 ./Qps;  dzs = 1 ./Qzs;

%---- Relationship between zero and pole freq.s -------
%---- GAMAs=Wps/Wzs and GAMAt=Wpt/Wzt -----------------

GAMAs =1.9;	%--- parallel filters with zeros aprox. 1 octave bellow poles
GAMAt = 0.6667 + (0.8667-0.6667)*xnorm;	 %--- cascade filters: poles bellow zeros
%GAMAt = 0.5+(0.85-0.5)*xnorm;

%--- Bilinear equivalent parameters: -----------------

Wzs = WCF./GAMAs;	%--- zero freq. for 2nd filters (fixed)
Wzt = WCTAP./GAMAt;	%--- zero freq. for cascade filters (fixed)

fizt = tan(Wzt/2);	%-- fixed
fizs = tan(Wzs/2);	%-- fixed
fipt = fiCTAP;		%-- for now! will be reestimated
fips = fiCF;

%------- filter coeficients (2nd order)---------------------

b0=fizt.^2; b1=dzt.*fizt;
a0=fipt.^2; a1=dpt.*fipt;
c0=1+a1+a0; d0=a0./(b0.*c0);

Bt=[ ((1+b1+b0).*d0)', (2*(b0-1).*d0)', ((1-b1+b0).*d0)' ];
At=[    ones(NS,1)   , (2*(a0-1)./c0)', ((1-a1+a0)./c0)' ];

b0=fizs.^2; b1=dzs.*fizs;
a0=fips.^2; a1=dps.*fips;
c0=1+a1+a0; d0=a0./(b0.*c0);

Bs=[ ((1+b1+b0).*d0)', (2*(b0-1).*d0)', ((1-b1+b0).*d0)' ];
As=[    ones(NS,1)   , (2*(a0-1)./c0)', ((1-a1+a0)./c0)' ];


%--- FIRST CASCADE FILTER --------------------------------------------------
%--- (different from others in order to have a broad bandwidth - to simulate 
%--- pseudoressonance). For this, the 1st cascade filter has only one zero:
%--- Model: bilinear transform of H(s)=(s*alpha*wp+wp^2)/(s^2+s*dp*wp+wp^2) 
%--- where wp/wz=alpha= const.
%--- dp=dpt(1); wp=wpt(1). Change only 1st cascade filter numerator

alpha=2.5/dpt(1);
b0=alpha*fipt(1)+fipt(1)^2;
b1=2*(fipt(1)^2);
b2=-alpha*fipt(1)+fipt(1)^2;
a0=1+dpt(1)*fipt(1)+fipt(1)^2;

Bt(1,:)=[b0, b1, b2]/a0;


%---- Middle Ear: -------------------------------------------
%---- 2nd order Butterworth High-Pass filter at 250Hz -------

%---- mudança para 350Hz

fi0=tan(pi*350/fs);
d=sqrt(2);
Bm=[1,-2,1];
Am=[1+d*fi0+fi0^2,2*(fi0^2-1),1-d*fi0+fi0^2];



%---- Earing Treshold ----------------------
%---- simulated in thresh.m

th = abs(Authresh(WCF*fs/2/pi)); %--- Fitting of the (inverted) earing treshold
Lth = th* 10^(45/20);	%--- adds +45dB
semilogx(fCF,dB(Lth))	%--- max|Yk(W)| should track this curve
hold on


%---- Middle ear response: --------

h=freqz(Bm,Am,W);



%---- PARAMETER REESTIMATION ---------------------------------------------
%---- force max|Pk(W)| to occurr at WCTAP chosing pole frequency fipt(k)
%---- force max|Yk(W)| to occurr at WCF with value (gain) Lth, by chosing 
%---- simultaneously dps(k) and fips(k)



for k=1:NS,		%--- for all sections...

   %----- find d/dW log(|P(k-1,W)|^2) and log(|P(k-1,W)|^2) ----------
   %----- at two frequencies of interest: WCF(k) and WCTAP(k) --------

   W2=[WCTAP(k);WCF(k)];		%--- (column vector)
   [dlogm,logm] = dlogmod2(Bm,Am,W2);
   
   for i=1:k-1,
     [dmi,mi] = dlogmod2(Bt(i,:),At(i,:),W2);
      dlogm = dlogm+dmi;
      logm  = logm+mi;			%--- update (log) deriv. and mod.
   end

   %--- Now, find fipt(k) such that |Pk(WCTAP)|=max(|Pk(W)|) -----

   [Bt(k,:),At(k,:)]=fip_max(Bt(k,:),At(k,:),dlogm(1),W2(1));


   %--- update log(|Pk(WCF)|^2) --------------------------

   [dmi,mi] = dlogmod2(Bt(k,:),At(k,:),W2);
   dlogm = dlogm+dmi;
   logm  = logm+mi;


   %--- find fips(k) and dps(k) such that |Yk(WCF)|=max(|Yk(W)|)=Lth -----
   %--- |Yk(Wc)|^2 = exp(logm + 2*log(|Hs|))=L^2 => |Hs| = L*exp(-logm/2)

   L=Lth(k)*exp(-logm(2)/2);
   As(k,:)=modA_Wc(Bs(k,:),dlogm(2),L,1,W2(2));
   a0=As(k,1);
   As(k,:)=As(k,:)/a0; Bs(k,:)=Bs(k,:)/a0;

   %--- plot responses ---------------------

   h=h.*freqz(Bt(k,:),At(k,:),W);
   y=h.*freqz(Bs(k,:),As(k,:),W);
   semilogx(f,db(y),f,db(h))
%%   semilogx([fCF(k),fCF(k)],[50,0],'r',[fCTAP(k),fCTAP(k)],[50,0],'b')

end

hold off


%--------------- SETS MAX GAIN 1 --------------


h=freqz(Bm,Am,W);

for k=1:NS,
   h=h.*freqz(Bt(k,:),At(k,:),W);
   hs=freqz(Bs(k,:),As(k,:),W);
   Y(k,:)=h.*hs;
end
%semilogx(f,db(Y))
maxY=max(max(abs(Y).'));
Bs=Bs/maxY;

%save coc52 Bm Am Bt At Bs As

