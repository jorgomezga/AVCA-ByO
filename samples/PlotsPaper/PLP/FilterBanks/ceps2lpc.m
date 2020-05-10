function [a,QME]=ceps2lpc(c)
% Convertion from lpc-ceptral coefficients to lpc coeficients.
%
% function [a,QME]=ceps2lpc(c)
%
% a	: colunm vector of lpc coeficients, with a(1)=1;
% QME	: lpc quadratic mean error 
% c	: ceptral coeficients c(1)..c(Q),   WITH  c(1)=log(QME)
%
% For matrices: 
%      c: cepstral coefficients in columns
%      a: lpc vectors (columns in a)
%
% Ref: J.Markel & A.Gray, "Linear Prediction of Speech", pp 230
% F. Perdigao, Coimbra, Aug. 96

QME=exp(c(1,:));
[Q,m]=size(c);
c(1,:)=[];
a=zeros(size(c));

a(1,:)=-c(1,:);
a(2,:)=-c(2,:)+(c(1,:).^2)/2;

for n=3:Q-1,
 k=(1:n-1);
 nkn=(n-k)'/n;
 a(n,:)=-c(n,:) - sum( a(k,:).*c(n-k,:).*nkn(:,ones(1,m)) ); 
end

%---- Note: the sum above must have at least 2 rows, i.e.: n>=3

a=[ones(1,m);a];