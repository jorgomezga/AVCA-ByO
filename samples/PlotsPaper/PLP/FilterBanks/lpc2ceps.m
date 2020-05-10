  function c=lpc2ceps(a,Q)
% function c=lpc2ceps(a,Q);
%
% Convertion from lpc coeficients to lpc-ceptral coeficients.
%
% a	: colunm vector of lpc coeficients, with a(1)=1;
% Q	: number of desired ceptral coeficients
% c	: ceptral coeficients c(1)..c(Q)
%
% For matrices: lpc vectors must be in columns in a
%
% c(n) = -a(n+1)-sum((n-k)*a(k+1)*c(n-k)/n), k=1..n-1 ,	for n<=p
% c(n) = -sum((n-k)*a(k+1)*c(n-k)/n), k=1..p , 		for n>p
%
% NOTE: c does not include the zeroth term: c(0)=log(QME)/2 
% where QME is the LPC prediction quadratic mean error.
%
% Ref: J.Markel & A.Gray, "Linear Prediction of Speech", pp 230

% F. Perdigao, Coimbra, Aug. 96
% fp@it.uc.pt



if Q<3 error('--> Q must be > 2'), end

[n,m]=size(a);

if Q<n, a=a(2:Q+1,:);
else    a=[a(2:n,:);zeros(Q-n+1,m)];	% zero-pad
end

c(1,:)=-a(1,:);
c(2,:)=-a(2,:)+(a(1,:).^2)/2;

for n=3:Q,
 k=(1:n-1);
 nkn=(n-k)'/n;
 c(n,:)=-a(n,:) - sum( a(k,:).*c(n-k,:).*nkn(:,ones(1,m)) ); 
end

%---- Note: the sum above must have at least 2 rows, i.e.: n>=3
