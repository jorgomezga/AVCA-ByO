function fi=fase(H,tol)
% Computs the unwraped phase of H(W) in multiples of pi.
% Calcula a fase de H(W) em multiplos de pi
%   fi=fase(H,tol)

if nargin==1,
 fi=unwrap(angle(H))/pi;
else
 fi=unwrap(angle(H),tol)/pi;
end