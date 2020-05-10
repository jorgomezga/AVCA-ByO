% Function to obtain the embedding dimension. First it tries to use the
% Cao's method, if it does not provide a value then the FNN method is used
% Inputs:
%       vSenal        Input signal
%       iTao          Delay time
%       iMaxDim       Maximum embedding dimension
%       iVerbosity    Verbosity level
%
% Outputs:
%       iDim:         Output embedding dimension

function iDim = CalcularEmbDim( vNormSignal, iTao, iMaxDim, iVerbosity )

if nargin<4
    iVerbosity=0;
end
if nargin<3
    iMaxDim=20; 
end
if nargin<2
    iTao= CalcularTao( vNormSignal );
end
if nargin<1
    error('Not enough input arguments')    
end

%%% Cao's method
iNsignal = length( vNormSignal );
if size( vNormSignal, 2 )>1
    vNormSignal = vNormSignal';
end
% TSTool constructor of type signal
oAux = signal( vNormSignal );   
% Cao method                               
vDim = cao_dim( oAux, iMaxDim-1, iTao, 4, iNsignal );   
% Convert from class signal to vector
vDimData = data( vDim );                                
vDimData = find( vDimData > 0.9 );
if ~isempty( vDimData )
    iDim = vDimData(1);
    if iVerbosity==1
        disp('                              ++Using Cao''s method')
    end
else
    %%% FNN if Cao's method do not converge
    if iVerbosity==1
        disp('                              ++Cao''s method did not converge, using FNN method')
    end
    de=1:length( vNormSignal );
    th=0.01;
    [iDim,~]=unfolding( vNormSignal, th, de, iTao );
end

if iDim < 2
    warning('                              ++Embedding dimension is less than 2, setting to 2.')
    iDim = 2;
end
if iDim > iMaxDim
    warning('                              ++iDim>iMaxDim. Setting tau to the value defined by iMaxTao')
    iDim = iMaxDim;
end