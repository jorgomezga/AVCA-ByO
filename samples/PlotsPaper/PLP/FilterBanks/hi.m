function h = Hi(ai,Qi,W)

u=0.5;                       % filtro de onda de pressao coclear
bb=0.5;                      % ref: Kates IEEE SP-39 No.12, p.2573 Dec.91
A=ai*ai*Qi;
B=ai*(1+u*Qi);
C=bb*u;

c=zeros(2,1);
c(1)=ai+u;
c(2)=ai-u;

%b=zeros(1,3);
%a=zeros(1,3);

b(1) = A+B+C;
b(2) = 2*(A-C);
b(3) = A-B+C;
a(1) = A+ai+Qi;
a(2) = 2*(A-Qi);
a(3) = A-ai+Qi;

bb=conv([ai ai],b);
aa=conv(c,a);
h = freqz(bb,aa,W);

