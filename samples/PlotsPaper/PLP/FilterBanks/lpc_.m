% LPC - Linear Prediction of order p
%
% - uses a optimized method for computing R when N=length(x) >> p
% 
% see lpc in signal toolbox
%
% function [a,EMQ,R]=lpc_(x,p)


function [a,EMQ,R]=lpc_(x,p)

x=x(:);
N=length(x);

ind=cumsum([(1:N);ones(p,N)]);	%--- indixes for computing autocorrelation R
R=zeros(size(ind));		%--- shape R

  y=[x;zeros(p,1)];		%--- x (column) plus p zeros
  R(:)=y(ind);			%--- reshape
  r=R*x;			%--- autocorrelation !

  a = [ 1; -toeplitz(r(1:p))\r(2:p+1)];	%--- LPC vector (column)
  EMQ=r'*a;			%--- quadratic mean error

