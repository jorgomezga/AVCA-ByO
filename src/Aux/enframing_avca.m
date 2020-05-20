function [mSignal]=enframing_avca( vSignal, iFrame, iSolape, sTipo)

%Function which performs framing and windowing. It permits to extract
%Glottal components and use them as input signal to be characterized
%
%   INPUTS:
%           - vSignal: audio signal vector
%           - iFs: sample frequency
%           - iFrame: frame length (in samples)
%           - iSolape: frame increment (in samples)
%           - sTipo: type of processing
%           - eOptions: structure containing special options. In this case,
%           if eOptions.removeSilence is true, silence frames are removed
%           using the vadsohn function from voicebox library.
%           - iVerbosity: if true, info about the process is printed on
%           screen.
%
%   OUPTUTS:
%           - mSignal: framed signal
%           - vOut:  Matrix containing the vocal tract model
%                                   for each frame
%
%

if nargin<4, error('Not enough input arguments!'); end


if any( sTipo=='R' )
    %Square window
    mSignal=enframe( vSignal, iFrame, iSolape );
elseif any ( sTipo=='N' )
    %Hanning window
    mSignal=enframe( vSignal, hanning( iFrame ), iSolape );
else %if any( sTipo=='M' ) || isempty( sTipo )
    %Hamming window
    mSignal=enframe( vSignal, hamming( iFrame ), iSolape );
end

end