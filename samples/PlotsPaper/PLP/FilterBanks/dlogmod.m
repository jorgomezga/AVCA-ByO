function y=dlogmod(ai,W,gama,Qi)
% computs d/dW ( ln |H(W)|^2 ) + gama
% as a function of ai in Kates' model
% see kates.m
%


while abs(W) > pi,
  W=W-2*pi;
end
if W<0,  W=-W;  end



u=0.5;
b=0.5;

ai2= ai^2;
u2 = u^2;
Qi2=Qi^2;
A = ai2*Qi;
B = ai*(1+u*Qi);
C = b*u;

c1 = ai2 + u2;
c2 = ai2 - u2;

rb0 = B*B + (A-C)^2;
rb1 = 2*(A*A - C*C);
rb2 = (A+C)^2 - B*B;

ra0 = ai2 + (A-Qi)^2;
ra1 = 2*(A*A - Qi2);
ra2 = (A+Qi)^2 - ai2;

sinW = sin(W);
cosW = cos(W);
sin2W = 2*sinW*cosW;

y = gama -sinW/(1+cosW) + ...
 c2*sinW/(c1+c2*cosW) -    ...
 (rb1*sinW + rb2*sin2W)/(rb0 + rb1*cosW + rb2*cosW*cosW) +  ...
 (ra1*sinW + ra2*sin2W)/(ra0 + ra1*cosW + ra2*cosW*cosW);

