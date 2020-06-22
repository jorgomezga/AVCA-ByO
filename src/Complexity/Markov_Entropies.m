function Entropies=Markov_Entropies(vSignal,iDim,iTao,rParam,flags)

% This function estimates a set of Markov Models-based entropies from a reconstructed 
% attractor using the embedding dimension iDim and the time delay iTao. 
% The estimator uses a discrete hidden Markov model whose codebook was set
% using the Subspace Constrained Mean Shift (SCMS) method for finding principal
% curves, which is used in this context as a "profile trajectory". 
%
% Please cite: Julián D. Arias-Londoño, Juan I. Godino-Llorente, "Entropies 
% from Markov Models as Complexity Measures of Embedded Attractors", 
% Entropy, vol. 17, no. 6, pp 3595-3620, 2015.
%
% Inputs:
%
%           vSignal: vector containing the time series to be analyzed.
%           iDim: Embedding dimension. Default 2.
%           iTao: Time delay. Default 1.
%           rParam: Window size for the kernel density estimator. Default 
%                   0.2*std(vSignal).
%			flags:  'M': Markov entropy Shannon and Renyi estimators
%				'm': Normalized Markov entropy Shannon and Renyi estimators				
%				'H': Conditional HMM entropy Shannon and Renyi estimators
%				'h': Normalized Conditional HMM entropy Shannon and Renyi estimators
%				'R': Recurrence state entropy Shannon and Renyi estimators
%				'r': Normalized Recurrence state entropy Shannon and Renyi estimators
%                   Defaul: 'MHR'
%
% Outputs: 
%
%           Entropies: Struct containing:
%
%           Set of entropies estimated over the attractor with dimension
%           iDim.
%
%           Entropies.OriginalDim.EMcs: Markov chain entropy:Shannon
%           Entropies.OriginalDim.EMcr: Markov chain entropy:Renyi 2
%           Entropies.OriginalDim.EhmmsN: Conditional Markov entropy (Normalized): Shannon
%           Entropies.OriginalDim.EhmmrN: Conditional Markov entropy
%                                         (Normalized): Renyi 2
%           Entropies.OriginalDim.EhmmsU: Conditional Markov entropy (Unnormalized): Shannon
%           Entropies.OriginalDim.EhmmrU: Conditional Markov entropy
%                                         (Unnormalized): Renyi 2
%           Entropies.OriginalDim.EMrds: Markov recurrence density: Shannon
%           Entropies.OriginalDim.EMrdr: Markov recurrence density: Renyi 2
%
%           ---------------------------------------------------------------
%
%           Set of entropies estimated as the average of the measures for 
%           the attractor with dimension iDim and iDim + 1.
%
%           Entropies.AverageDim.EMcs: Markov chain entropy:Shannon
%           Entropies.AverageDim.EMcr: Markov chain entropy:Renyi 2
%           Entropies.AverageDim.EhmmsN: Conditional Markov entropy (Normalized): Shannon
%           Entropies.AverageDim.EhmmrN: Conditional Markov entropy
%                                         (Normalized): Renyi 2
%           Entropies.AverageDim.EhmmsU: Conditional Markov entropy (Unnormalized): Shannon
%           Entropies.AverageDim.EhmmrU: Conditional Markov entropy
%                                         (Unnormalized): Renyi 2
%           Entropies.AverageDim.EMrds: Markov recurrence density: Shannon
%           Entropies.AverageDim.EMrdr: Markov recurrence density: Renyi 2


if nargin < 1, error('Incorect number of parameter'); end
if nargin < 2, iDim = 2; end
if nargin < 3, iTao = 1; end
%------------------- Preproceso -------------------------------------------
if size(vSignal,2) > 1
    vSignal = vSignal';
end
vSenal2 = vSignal-mean(vSignal);
vSenal2 = (vSenal2 - min(vSenal2))./(max(vSenal2)+abs(min(vSenal2))+ 0.01);
%--------------------------------------------------------------------------
if nargin < 4, rParam = 0.2*std(vSenal2); end
if nargin < 5, flags = 'MHR'; end
%--------------------------------------------------------------------------
[mAtractor]=embeb(vSenal2,iDim,iTao);
[Salida,~]=Conditional_DHMM_KernelCorrEntropy2(mAtractor,size(mAtractor,1),2,rParam,flags);
[mAtractor]=embeb(vSenal2,iDim + 1,iTao);
[Salida2,~]=Conditional_DHMM_KernelCorrEntropy2(mAtractor,size(mAtractor,1),2,rParam,flags);

if any( flags=='M')
    Entropies.OriginalDim.EMcs = Salida.Emcs;
    Entropies.OriginalDim.EMcr = Salida.Emcr;
end
if any( flags=='H' )
    Entropies.OriginalDim.EhmmsN = Salida.Ehmms;
    Entropies.OriginalDim.EhmmrN = Salida.Ehmmr;
    Entropies.OriginalDim.EhmmsU = Salida.Ehmms2;
    Entropies.OriginalDim.EhmmrU = Salida.Ehmmr2;
end
if any( flags=='R' )
    Entropies.OriginalDim.EMrds = Salida.MRDs;
    Entropies.OriginalDim.EMrdr = Salida.MRDr;
end
%------------------------------------------------------------------

%------------------------------------------------------------------
if any( flags=='m' )
    Entropies.AverageDim.EMcs = (Salida2.Emcs + Salida.Emcs )/2;
    Entropies.AverageDim.EMcr = (Salida2.Emcr + Salida.Emcr )/2;
end
if any( flags=='h' )
    Entropies.AverageDim.EhmmsN = (Salida.Ehmms + Salida2.Ehmms)/2;
    Entropies.AverageDim.EhmmrN = (Salida.Ehmmr + Salida2.Ehmmr)/2;
    Entropies.AverageDim.EhmmsU = (Salida.Ehmms2 + Salida2.Ehmms2)/2;%Sin normalizar
    Entropies.AverageDim.EhmmrU = (Salida.Ehmmr2 + Salida2.Ehmmr2)/2;%Sin normalizar
end
if any( flags=='r' )
    Entropies.AverageDim.EMrds = (Salida2.MRDs + Salida.MRDs)/2;
    Entropies.AverageDim.EMrdr = (Salida2.MRDr + Salida.MRDr)/2;
end
