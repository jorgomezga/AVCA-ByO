function[mAtractor]=embeb(vSignal,iDim,iTao)

% Reconstruct state space vectors in a iDim dimensional space, using iTao delays
% Input:
%		vSignal:	Input signal
%		iDim:		Embedding dimension
%		iTao:		delay
%
% Output:
%		mAtractor: 	State space reconstruction. Matrix of size (iDim) x (N-(iDim)*iTao)


mAtractor=vSignal((1:length(vSignal)-((iDim-1)*iTao))'*ones(1,iDim)+ones(length(vSignal)-((iDim-1)*iTao),1)*(0:iDim-1)*iTao);