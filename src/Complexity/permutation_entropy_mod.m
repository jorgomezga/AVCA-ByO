function [PE, PEn] = permutation_entropy_mod( X )

% [PE, PEn] = permutation_entropy(y, m, t)
% 
% This function calculates the permutation entropy of the input column
% vector.
% 
% X     - Atractor
% PE    - permutation entropy
% PEn   - normalized permutation entropy
% 
% Implemented according to:
% Zanin, M.; Zunino, L.; Rosso, O.A.; Papo, D. Permutation Entropy and Its
% Main Biomedical and Econophysics Applications: A Review. Entropy 2012,
% 14, 1553-1577.

%% Get the patterns
[~, pat] = sort(X.');

%% Get the probabilities
[~, ~, n] = unique(pat.','rows');

m = size( X, 2 );
m_fac = factorial(m);

p = hist(n,max(n));
p = p./length(n);

% %% Calculate the permutation entropies
PE = -sum(p.*log2(p));
PEn = PE/log2(m_fac);
