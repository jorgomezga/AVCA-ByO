function SndF = SecondF(Wc,Qp,Q0,W)

betap = tan(Wc/2);
beta0 = betap/2;
b02 = beta0*beta0;
bp2 = betap*betap;
b0Q = beta0/Q0;
bpQ = betap/Qp;

b(1) = 4*(b02 + b0Q + 1);
b(2) = 8*(b02-1);
b(3) = 4*(b02 - b0Q + 1);

a(1) = bp2 + bpQ + 1;
a(2) = 2*(bp2-1);
a(3) = bp2 - bpQ +1;

SndF = freqz(b,a,W);
