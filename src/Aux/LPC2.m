function [vA,rG] = LPC2( vFrame, iOrderP, iPlot )

% Computes the coeffiecients vA and the gain of the filter using a LPC
% analysis of a voice signal
%       A = {a0=1, -a1, -a2,..., -ap}
%
% Inputs:
%   vFrame               = Input vFrame
%   iOrderP              = Prediction order
%   iPlot                = logical value. If true, the coefficients are
%                           plotted
% Outputs:
%   vA        = Coefficients of the LPC analysis
%   rG        = Gain of the filter in the LPC analysis

if iscolumn( vFrame )
    vFrame = vFrame';
end

vK = zeros(1,iOrderP);
va = zeros(iOrderP,iOrderP);
vE = zeros(1,iOrderP);

% Calculation of the normalized autocorrelation
vR = xcorr(vFrame,'biased'); % biased es para normalizar a la unidad
iM = length(vFrame);
vR(1:iM-1) = [];

% Durbin iterative method. To avoid troubles in the calculations, the first
% iteration is made outside the loop
e=vR(1);
vK(1)=-(vR(2)/e);
va(1,1)=vK(1);
vE(1)=(1-vK(1)^2)*e;

% MARKHOUL implementation (38a)->(38e)
for i=2:iOrderP
    vK(i)=-( ( vR(i+1) + va((1:(i-1)),(i-1))' * (fliplr(vR(2:(i-1+1))))' ) / vE(i-1) ) ;
    va(i,i)=vK(i);
    for j=1:(i-1)
        va(j,i)=va(j,i-1) + vK(i)*va(i-j,i-1);
    end
    vE(i)=(1-vK(i)^2)*vE(i-1);
end

% a0 is always 1, and the remaining elements (-a1, -a2, ..., -ap)  are the
% last iterations (last column).
vA=[1 va((1:iOrderP), iOrderP)'];
% rG^2=Ep since if the input signal is an impulse function, the error will
% be the energy of the impulse
rG = sqrt(vE(iOrderP));

% The advantage of this method in addition to its speed (although Matlab uses
% the same method) is that it does not only provides the coefficients
% of the prediction of order p but also those of order p-1, p-2, ...
% It also gives us the Error E = SUM (en ^ 2), being en = sn-sn~
% We can see how by increasing the order of the prediction
% the approximation is gets better and the Error gests smaller.
if iPlot
    figure(1);
    iSizeWin=256;
    
    iOrderP=3;
    a_p=[1 va((1:iOrderP), iOrderP)'];
    [H,~]=freqz(iSizeWin.^.5*1,a_p,iSizeWin/2);
    env=log(abs(H));
    subplot(2,2,1), plot(env);
    
    iOrderP=6;
    a_p=[1 va((1:iOrderP), iOrderP)'];
    [H,~]=freqz(iSizeWin.^.5*1,a_p,iSizeWin/2);
    env=log(abs(H));
    subplot(2,2,2), plot(env);
    
    iOrderP=12;
    a_p=[1 va((1:iOrderP), iOrderP)'];
    [H,~]=freqz(iSizeWin.^.5*1,a_p,iSizeWin/2);
    env=log(abs(H));
    subplot(2,2,3), plot(env);
    
    iOrderP=13;
    a_p=[1 va((1:iOrderP), iOrderP)'];
    [H,~]=freqz(iSizeWin.^.5*1,a_p,iSizeWin/2);
    env=log(abs(H));
    subplot(2,2,4), plot(env);
    
    iOrderP=49;
    a_p=[1 va((1:iOrderP), iOrderP)'];
    [H,~]=freqz(iSizeWin.^.5*1,a_p,iSizeWin/2);
    env=log(abs(H));
    subplot(2,2,4), plot(env);
    
    figure(2);
    plot(vE(2:length(vE)));
    title('Error signal');
    
end