function rEp=ErrorPredNorm( vFrame, iOrderLPC )

% Normalized Prediction error in dB
%
% Inputs:
%   vFrame               = Input vFrame
%   iOrderLPC            = Prediction order
% Outputs:
%   rEp                  = Normalized prediction error

iLengthFrame = length( vFrame );
if length( vFrame )-iOrderLPC < iOrderLPC+1 
    error('The length is so short or the order of the LPC analysis is too high'); 
end

% LPC computation: vA=[1, -a1, -a2, ..., -ap]
[vA,~,~] = lpcauto( vFrame(iOrderLPC+1:end), iOrderLPC );

% Log Energy of the frame
rEs = LogEnergy( vFrame(iOrderLPC+1:end) );

vPhi = zeros(1, iOrderLPC-1);
for k=0:iOrderLPC
   rSuma=0;
   for j=iOrderLPC+1:iLengthFrame
      rSuma=rSuma+vFrame(j)*vFrame(j-k);
   end
   vPhi(k+1)=rSuma/(iLengthFrame-iOrderLPC);
end

rSuma=0;
for i=2:iOrderLPC+1
   rSuma=rSuma+vA(i)*vPhi(i);
end

rEp = rEs-10*log10( eps+abs(rSuma+vPhi(1)) );