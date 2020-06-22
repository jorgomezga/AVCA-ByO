% Modificada sobre la funcion de Max Little
% Performs fast recurrence period density entropy (RPDE) analysis on an input signal to
% obtain an estimate of the H_norm value.
%
% Useage:
% [H_norm, rpd] = rpde(x)
% [H_norm, rpd] = rpde(x, epsilon)
% [H_norm, rpd] = rpde(x, epsilon, T_max)
% Inputs
%    x       - Embedding Matrix (Must be in format: n rows x d dimensions)
% Optional inputs
%    epsilon - recurrence neighbourhood radius
%              (If not specified, then a suitable value is chosen automatically)
%    T_max   - maximum recurrence time
%              (If not specified, then all recurrence times are returned)
% Outputs:
%    H_norm  - Estimated RPDE value
%    rpd     - Estimated recurrence period density
%
% (c) 2007 Max Little. If you use this code, please cite:
% Exploiting Nonlinear Recurrence and Fractal Scaling Properties for Voice Disorder Detection
% M. Little, P. McSharry, S. Roberts, D. Costello, I. Moroz (2007),
% BioMedical Engineering OnLine 2007, 6:23

function [H_norm, rpd] = rpde_matlab_C(embedSequence, epsilon, T_max)

if ((nargin < 1) || (nargin > 3))
    help rpde;
    return;
end

if (nargin < 2)
    epsilon = 0.12;
end

if (nargin < 3)
    T_max = -1;
end

embedDim=size(embedSequence,2);             %Numero de dimensiones
embedSequence=reshape(embedSequence',[],1); %Se concatenan las filas de la matriz
                                            %Esto para que quede en el
                                            %formato que acepta la funcion
                                            %de close_ret_mod

%Funcion que halla los close returns de los vectores de estado                                           
rpd = close_ret_mod(embedSequence, embedDim ,epsilon);

if (T_max > -1)
    rpd = rpd(1:T_max);
end
rpd = rpd/( sum(rpd) + eps );

N = length(rpd);
H = 0;
for j = 1:N
   H = H - rpd(j) * logz(rpd(j));
end
H_norm = H/log(N);


function y = logz(x)
if (x > 0)
   y = log(x);
else
   y = 0;
end
