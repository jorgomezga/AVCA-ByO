function [AE, SE, mSE, GE, FE] = CalculateRegularity( mAttractorp1, mAttractor,...
    rParam, sTypeNorm, bFast )

% Calculates a series of regularity features from the inputted attractor
% mAttractorp1 of dimension (m+1), using the tolerance value rParam.
%
% Approximate entropy (AE): Pincus, S. M., Gladstone, I. M., & Ehrenkranz, R. A. (1991).
%       A regularity statistic for medical data analysis.
%       Journal of clinical monitoring, 7(4), 335-345.
% Sample entropy (SE): Richman, J. S., & Moorman, J. R. (2000).
%       Physiological time-series analysis using approximate entropy and sample entropy.
%       American Journal of Physiology-Heart and Circulatory Physiology, 278(6), H2039-H2049.
% modified Sample Entropy (mSE): Xie, H. B., He, W. X., & Liu, H. (2008).
%       Measuring time series regularity using nonlinear similarity-based sample entropy.
%       Physics Letters A, 372(48), 7140-7146.
% Gaussian Sample Entropy (GE): Xu, L. S., Wang, K. Q., & Wang, L. (2005, August).
%       Gaussian kernel approximate entropy algorithm for analyzing irregularity of time-series.
%       In 2005 international conference on machine learning and cybernetics (Vol. 9, pp. 5605-5608). IEEE.
% Fuzzy Entropy (FE): Chen, W., Zhuang, J., Yu, W., & Wang, Z. (2009).
%       Measuring complexity using fuzzyen, apen, and sampen.
%       Medical engineering & physics, 31(1), 61-68.

% Inputs:
%   mAttractorp1 = Embedded matrix of dimension m+1. Size n x (m+1)
%   mAttractor   = Embedded matrix of dimension m. Size n x m
%   rParam       = Tolerance parameter for the calculation of the distances
%   sTypeNorm    = Type of norm in the calculation of distances. View pdist
%                  function in Matlab for options. Typical values include:
%                  'chebychev': inf-norm
%                  'euclidean': 2-norm
%                  'cityblock': 1-norm
%   bFast        = boolean indicating if using the fast or the slow
%                  version. The fast version uses pdist, hence requires a
%                  larger amount of ram. The slow version does not have
%                  large ram requirement
% Outputs:
%   AE  = Approximate Entropy
%   SE  = Sample Entropy
%   mSE = modified sample Entropy
%   GE  = Gaussian sample Entropy
%   FE  = Fuzzy Entropy

if nargin < 5, bFast = true;            end
if nargin < 4, sTypeNorm = 'chebychev'; end
if nargin < 3, rParam = 0.1;            end
if nargin < 2, mAttractor  = mAttractorp1(:,1:end-1);   end
if nargin < 1, error( 'Not enough input parameters!' ); end

if bFast
    % Fast implementation -> Using pdist, faster but requiring more RAM
    [AE, SE, mSE, GE, FE] = FastRegularity( mAttractorp1, mAttractor, rParam, sTypeNorm );
else
    % Slow implementation -> Using loops, slower but not requiring as much RAM as
    % with pdist
    [AE, SE, mSE, GE, FE] = SlowRegularity( mAttractorp1, mAttractor, rParam, sTypeNorm );
end

end


%% Fast implementation using pdist
function [AE, SE, mSE, GE, FE] = FastRegularity( mAttractorp1, mAttractor, rParam, sTypeNorm )

% Distance calculation
% Size of the embedding matrix
[nfil,~]   = size( mAttractor );
[nfilp1,~] = size( mAttractorp1 );

% Dim m
dist_p  = pdist( mAttractor, sTypeNorm );
% Dim m+1
dist_p1 = pdist( mAttractorp1, sTypeNorm );

% ApEn
AE     = ApEn_Aux( mAttractorp1, mAttractor, rParam, sTypeNorm );
% SampEn
SEm    = sum( SampEn( dist_p, rParam )/nfil )/nfil ;
SEmp1  = sum( SampEn( dist_p1, rParam )/nfilp1 )/nfilp1;
SE     = -log( SEmp1 / SEm );
% mSampEn
mSEm   = sum( mSampEn( dist_p, rParam )/nfil )/nfil;
mSEmp1 = sum( mSampEn( dist_p1, rParam )/nfilp1 )/nfilp1;
mSE    = -log( mSEmp1 / mSEm );
% GSampEn
GEm    = sum( GSampEn( dist_p, rParam )/nfil )/nfil;
GEmp1  = sum( GSampEn( dist_p1, rParam )/nfilp1 )/nfilp1;
GE     = -log( GEmp1 / GEm );
% FuzzyEn
re = 2;
FEm    = sum( FuzzyEn( dist_p, rParam, re )/nfil )/nfil; 
FEmp1  = sum( FuzzyEn( dist_p1, rParam, re )/nfilp1 )/nfilp1;
FE     = -log( FEmp1 / FEm );

end


%% Slow implementation using loops
function [AE, SE, mSE, GE, FE] = SlowRegularity( mAttractorp1, mAttractor, rParam, sTypeNorm )

re = 2;
% Size of the embedding matrix
[nfil,~]   = size( mAttractor );
[nfilp1,~] = size( mAttractorp1 );

% ApEn
AEm    = zeros( 1, nfil );
AEmp1  = zeros( 1, nfilp1 );
% SampEn
SEm    = zeros( 1, nfil );
SEmp1  = zeros( 1, nfilp1 );
% mSampEn
mSEm   = zeros( 1, nfil );
mSEmp1 = zeros( 1, nfilp1 );
% GSampEn
GEm    = zeros( 1, nfil );
GEmp1  = zeros( 1, nfilp1 );
% FuzzyEn
FEm    = zeros( 1, nfil );
FEmp1  = zeros( 1, nfilp1 );

for i=1:nfil    
    % Dimension m:    
    dist_p = pdist2( mAttractor, mAttractor(i,:), sTypeNorm );
        
    % ApEn
    AEm(i)    = log( sum( SampEn( dist_p, rParam) ) /nfil );
    % SampEn
    SEm(i)    = sum( SampEn( dist_p(i+1:nfil,:), rParam) )/nfil;
    % mSampEn
    mSEm(i)   = sum( mSampEn( dist_p(i+1:nfil,:), rParam ) )/nfil;
    % GSampEn
    GEm(i)    = sum( GSampEn( dist_p(i+1:nfil,:), rParam ) )/nfil;
    % FuzzyEn
    FEm(i)    = sum( FuzzyEn( dist_p(i+1:nfil,:), rParam, re) )/nfil;    
    
    % Dimension m+1
    if i<=nfilp1        
        dist_p1 = pdist2( mAttractorp1, mAttractorp1(i,:), sTypeNorm );
        
        % ApEn
        AEmp1(i)  = log( sum( SampEn( dist_p1, rParam ) ) /nfilp1 );
        % SampEn
        SEmp1(i)  = sum( SampEn(  dist_p1(i+1:nfilp1,:), rParam ) )/nfilp1;
        % mSampEn
        mSEmp1(i) = sum( mSampEn( dist_p1(i+1:nfilp1,:), rParam ) )/nfilp1;
        % GSampEn
        GEmp1(i)  = sum( GSampEn( dist_p1(i+1:nfilp1,:), rParam ) )/nfilp1;
        % FuzzyEn
        FEmp1(i)  = sum( FuzzyEn( dist_p1(i+1:nfilp1,:), rParam, re ) )/nfilp1;        
    end
end

% ApEn
AEm   = mean( AEm );
AEmp1 = mean( AEmp1 );
if AEm==0
    AEm=eps;
elseif AEmp1==0
    AEmp1=eps;
end
AE    = AEm - AEmp1;
% SampEn
SEm   = mean( SEm );
SEmp1 = mean( SEmp1 );
SE    = -log( SEmp1 / SEm );
% mSampEn
mSEm   = mean( mSEm );
mSEmp1 = mean( mSEmp1 );
mSE    = -log( mSEmp1 / mSEm );
% GSampEn
GEm   = mean( GEm );
GEmp1 = mean( GEmp1 );
GE    = -log( GEmp1 / GEm );
% FuzzyEn
FEm   = mean( FEm );
FEmp1 = mean( FEmp1 );
FE    = -log( FEmp1 / FEm );

end


%% Membership functions
function count = SampEn( dist, rn )
    count = (dist<= rn);
end

function count = GSampEn( dist, rn )
    count = exp(-(dist.^2)./(10*rn^2));
end

function count = mSampEn( dist, rn )
    count = 1./(1+exp(abs((dist+0.5)/rn)));
end

function count = FuzzyEn( dist, rn, re )
    count = exp( -(dist./rn).^re );
end


%% Auxiliary function to calculate ApEn
function AE = ApEn_Aux( mAttractorp1, mAttractor, rParam, sTypeNorm )

% Dimension m
[nfil,~]   = size( mAttractor );
AEp        = zeros( 1, nfil );
% Dimension m+1
[nfilp1,~] = size( mAttractorp1 );
AEmp1      = zeros( 1, nfilp1 );

for i=1:nfil
    % Dimension m
    dist_m = pdist2( mAttractor, mAttractor(i,:), sTypeNorm );
    AEp(i) = log( sum( feval( 'SampEn',  dist_m, rParam )/nfil ) );
    
    % Dimension m+1
    if i<=nfilp1
        dist_m1  = pdist2( mAttractorp1, mAttractorp1(i,:), sTypeNorm );
        AEmp1(i) = log( sum( SampEn( dist_m1, rParam )/nfilp1 ) );
    end
end

AEp   = mean(AEp);
AEmp1 = mean(AEmp1);
if AEp==0
    AEp=eps;
elseif AEmp1==0
    AEmp1=eps;
end
AE= AEp-AEmp1;

end