 function vParameters = Fluctuation( vSignal, iFs )

% Calculation of the perturbation and fluctuation parameters
%
% Inputs:
%   vSignal              = Input signal
%   iFs                  = Sampling frequency
%
% Outputs:
%   vParameters: Vector contatenated with the following features:
%           rFftr:      frequency of the low-frequency component of higher intensity that modulates
%                       the pitch period.
%           rFatr:      frequency of the low-frequency component of higher intensity that modulates
%                       the pitch amplitude.
%           rFTRI:      (Fo Tremor Intensity Index)
%           rATRI:      (Amplitude Tremor Intensity Index)

%% Fluctuation
try
    [rFTRI, rATRI, rFftr, rFatr]=Tremor( vSignal, iFs );
catch
    warning('Tremor parameters not calculated');
   rFTRI=NaN ;
   rATRI=NaN ;
   rFftr=NaN ;
   rFatr=NaN ;
end

vParameters = [rFTRI, rATRI, rFftr, rFatr];