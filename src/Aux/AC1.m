function rC1=AC1( vFrame )

% Normalized autocorrelation coefficient for a lag of 1 of a vFrame
%
% Inputs:
%   vFrame               = Input vFrame
% Outputs:
%   rC1                  = Normalized autocorrelation coefficient for a lag of 1

rNumerator=0;
for n=2:length( vFrame )
    rNumerator = rNumerator+vFrame(n)*vFrame(n-1);
end

rProd1 = sum(vFrame( 2:length( vFrame ) ).^2);
rProd2 = sum(vFrame( 1:length( vFrame )-1 ).^2);
rDenominator = sqrt( rProd1*rProd2 );

rC1 = rNumerator/rDenominator;