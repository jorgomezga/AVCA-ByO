% Triangular function centered at T0 and width (Tupp-Tlow)
% given t
%         
%       /\1
%      /  \
%     /    \
%    /      \
%   /        \
%__/__________\_____
% Tlow  T0    Tupp   t
%
% function f = triangf(t,T0,Tlow,Tupp)

function f = triangf(t,T0,Tlow,Tupp)


f=(t-Tlow)/(T0-Tlow).*(t>=Tlow & t<=T0) + ...
  (Tupp-t)/(Tupp-T0).*(t>T0 & t<=Tupp);

