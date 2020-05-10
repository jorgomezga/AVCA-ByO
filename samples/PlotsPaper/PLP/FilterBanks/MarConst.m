% function [A,B,Vcenter]=MarConst(Th,fspo,fsat,n)
%
% Calcula as constantes A,B do modelo de Martens-Immerseel
% de forma a obter um modelo de uma fibra de limiar
% de excitação Th (em dB) e com taxas de disparos espontâneas e 
% de saturação, fspo e fsat, respectivamente.
% As outras constantes do modelo não são alteradas de forma
% a não alterar significativamente as constantes de tempo
% de adaptação.
% Devolve também Vcenter, a amplitude da sinosoide (em dB) que 
% conduz a uma taxa de disparos f=(fsat+fspo)/2.
%

function [A,B,Vcenter]=MarConst(Th,fspo,fsat,n)

A=10^(Th/20);
B = (A^(1/n))*( (fsat/fspo)^(1/n) - 1 );



%--- RATE-INTENSITY curve --------

Vdb= linspace(Th-20,Th+80,200);
V=10.0.^(Vdb/20);
i=find(V<A);
uv(i)=A*ones(size(i));
i=find(V>=A);
uv(i)=(A*acos(-A./V(i))+sqrt(V(i).^2-A^2))/pi;
uf=fsat*uv./((B+uv.^(1/n)).^n);
plot(Vdb,uf,'c')

%------- Dynamic range in 5% to 95% of firing rate

v1=zcross(uf-1.05*fspo,Vdb);
v2=zcross(uf-0.95*fsat,Vdb);

disp('Martens-Immerseel IHC/Synapse model.');
disp(['A=',num2str(A,3),' B=',num2str(B,3)]);
disp(['Dynamic range (5% to 95% of absolute values): ', num2str(v2-v1,3), ' dB']);
disp(['Rate threshold: ', num2str(v1,3), ' dB']);
disp(['Sat. threshold: ', num2str(v2,3), ' dB']);


fcenter = (fsat+fspo)/2;
%i=find_ind(uf,fcenter);
%Vcenter=Vdb(i);

v=( B/((fsat/fcenter)^(1/n) -1) )^n;
Vcenter=db(sqrt(A^2+(pi*(v-A/2))^2))

