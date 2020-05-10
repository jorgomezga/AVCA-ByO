  function D=deltas(X,N)
% function D=deltas(X,N)
%
% Deltas: linear regression of vectors (e.g. cepstral vectors)
% X : matrix with columns as feature vectors: X(Nbank,Nframes)
% N : window size: uses values from X(:,n-N) to X(:,n+N) (2N+1) values
%     (default: N=2)
% D : delta matrix (with the same size as X)
%
% D(k)= G*sum(n*X(n+k)), n=-N:N
% G=1/sqrt(sum(n^2)), n=-N:N - (normalizing constant)
%
% NOTE: the first and last N values are copies from the adjacent calculated deltas. 

% F. Perdigao, IT-Coimbra, Aug. 96 (fp@it.uc.pt)


if nargin==1, N=2; G=0.3162278; 
else
 G= 1/sqrt(2*sum((1:N).^2));	%--- or: G=1/sqrt(N*(N+1)*(2*N+1)/3);
end

[Q,Nf]=size(X);
if Nf < (2*N+1), error('too few columns!'), end

N1=N+1;		%--- 1st vector index for full evaluation
Nf1=Nf-N;	%--- last vector index  "   "     "
D=0;

for n=1:N,
 D=D+n*( X(:,N1+n:Nf1+n) - X(:,N1-n:Nf1-n) );
end

Di=D(:,1);		%--- puts N samples in the begin and the end
Df=D(:,Nf-2*N);		%--- in order to obtain Nf samples as in X

D=G*[Di(:,ones(1,N)),D,Df(:,ones(1,N))];


