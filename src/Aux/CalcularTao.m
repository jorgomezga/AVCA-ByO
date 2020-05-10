% Function to obtain the delay time (tau). First it tries to use the automutual
% information criterion, if it does not provide a value then the first zero
% of the autocorrelation function is used
% Inputs:
%       vSenal        Input signal
%       iMaxTao       Maximum delay time
%       iVerbosity    Verbosity level
%
% Outputs:
%       iTao:         Output delay time

function iTao = CalcularTao( vSenal, iMaxTao, iVerbosity )

if nargin<3
    iVerbosity=0;
end
if nargin<2  
    iMinDim = 2;
    iMaxTao = floor( length( vSenal )/iMinDim ); 
end
if nargin<1
    error('Not enough input arguments')    
end

iNsignal = length( vSenal );
if size( vSenal, 2 )>1
    vSenal = vSenal';
end
    
% Grouping intervales, HON page 198. 
k = floor( 5*log10( iNsignal ) );   
% Automutual information (TSTOOL)
vTao = amutual( vSenal, iMaxTao, k );          
vAux = diff( vTao );
vAux = find( vAux > 0 );
if ~isempty( vAux )
    if iVerbosity==1
        disp('                              ++Tau using the automutual information criterion')
    end
    iTao = vAux(1)-1;                      
else
    %%% First zero of the autocorrelation as the amutual criterion did not
    %%% converge
    if iVerbosity==1
        disp('                              ++Automutual information criterion did not converge, use the first zero of the autocorrelation function')
    end
    iTao = firstzero( vSenal );
end    

if iTao>iMaxTao
    iTao = iMaxTao;
    warning('                              ++iTao>iMaxTao. Setting tau to the value defined by iMaxTao')
end