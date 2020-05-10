% kates
% Computs Kates' Auditoy Model freq. responses
% Ref: J. Kates, "A Time-Domain Digital Cochlear Model", 
% IEEE, Signal Proc., Vol.39, No.12, pp. 2573-2592, Dec. 1991
% Note: needs Signal Processing Toolbox (freqz.m)


clear
clc
disp('Kates´ model Calculation...')
disp('please wait (about 1 min)...')

fs=40000;
NS=112;
NPOINTS=512;                      % no. points in freq.
% --------------------

x  = linspace(0.9540228,0.0732451,NS);	%--- CFs from 16kHz to 100Hz
fc = 160*(10 .^ (2.1*x) -0.8);    	% Liberman map
Wc = 2*pi*fc/fs;

f  = logspace(2,4,NPOINTS);
W  = 2*pi*f'/fs ;                 % for freq. responses.

Qmax = 0.217*x + 0.264;           % (pp.2576 1st column)
Qp = 1.5*(1+fc/1000);             % Quality factors for 2nd filter
Q0 = 2*Qp;

G=zeros(NPOINTS,NS);
V=G;
U=V;
%----------------------------------------------------------------
% 1st section

   k=1;
   Wk = Wc(1);
   Q=Qmax(1);
                                  % guess = Wc
   a(1) =  fzero2('dlogmod',Wk,Wk,0,Q);
   G(:,1) = Hi(a(1),Qmax(1),W);   % G(k=1)

k0   = tan(a(1)/4);                   % ????
HP1  = freqz([1,-1],[k0+1,k0-1],W);   % HPF 1st order
SndF = SecondF(Wk,Qp(1),Q0(1),W);

V(:,1) = G(:,1) .* HP1;               % velocidade
U(:,1) = V(:,1) .* SndF;              % 2nd filter output

%semilogx(f,dB(U(:,1)), f, dB(G(:,1)))
%semilogx(f,dB(U(:,1)))
%hold on



%------------ CYCLE --------------------------------------

    for k=2:NS,

        gama = 0;
        Wk = Wc(k);

            for i=1:(k-1),
              gama = dlogmod(a(i),Wk,gama,Qmax(i));
            end

        guess = Wk;
        a(k) = fzero2('dlogmod',guess,Wk,gama,Qmax(k));
        H = Hi(a(k),Qmax(k),W);
        G(:,k) = H .* G(:,k-1);

        k0=tan(a(k)/4);
        HP1  = freqz([1,-1],[k0+1,k0-1],W);
        SndF = SecondF(Wk,Qp(k),Q0(k),W);

        V(:,k) = G(:,k) .* HP1;                 % velocity
        U(:,k) = V(:,k) .* SndF;                % 2nd filter output

    if (k == 30)|(k==50)|(k==70)|(k==90)|(k==100),
        %semilogx(f,20*log10(abs(U(:,k))), f, dB(G(:,k)))
        semilogx(f,db(U(:,k),1e-3)), hold on
        drawnow
    end
 end

hold off
axis([1e2,1e4,-60,60])


disp(' Reprodution of Kates´ paper figure, magnitude')

semilogx(f,db(U(:,[30,50,70,90,100]),1e-3))
axis([1e2,1e4,-60,60])
pause


disp(' Phase ')

for k=[30,50,70,90,100],
 z=fase(U(:,k));
 i=find(z<-8);
 z(i)=-8*ones(size(i));
 semilogx(f,z), hold on
end
hold off

disp(' Now variable "U" have 112 tap filter responses, each with 512 points in freq.') 

