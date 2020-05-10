close all
clc

t = -1:1/1000:1;
r = 0.2;

%% GapEn
GapEn = exp( -t.^2/(10*r^2) );
figure
plot( t, GapEn )
axis off

%% Fuzzy Entropy
Fuzzy = exp( -t.^2/r );
figure
plot( t, Fuzzy )
axis off

%% mSampEn
mSampEn1 = 1./ (1+exp( (t-0.5)/r ));
mSampEn2 = 1./ (1+exp( (-t-0.5)/r ));
figure
plot( [mSampEn2,mSampEn1] )
axis off

%% Heaviside
tn1 = zeros( 1, 300 );
tn2 = ones( 1, 600 );

HS = [tn1, tn2, tn1];
plot( HS )
axis off