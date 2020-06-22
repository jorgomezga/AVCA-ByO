function [P, i_P, mPicosPos]=Picos( vSignal, iFs, iInicio, iFinal, vPitchPeriods, pos_neg )

% Computes the amplitudes P and the i_P positions of the positive or negative peaks
% (depending on the value of the pos_neg parameter), cycle by cycle, of a voice frame.
% It also returns the peaks, with the same duration as vSignal, which contains
% the peaks found in the analysis.
%
% Input parameters
%   vSignal:        is a column vector containing the complete speech signal
%   iFs:            is the sampling frequency
% 	iInicio:        is the first sample of the voice to be analyzed
%   iFinal:         is the last sample of the section to be analyzed
%   vPitchPeriods:  is the sequence of pitch period values (in samples)
%                   calculated on windows of duration 3To without overlap,
%                   where To is the pitch period in samples from the previous segment,
%                   and 34ms the duration of the first window, as Boyanov does.
%                   It is the result of applying the pitch function.
%   pos_neg:        indicates which peaks are going to be obtained:
%                   pos_neg = 1: positive
%                   pos_neg = -1: negative
%
% Output parameters
%
% P:                The amplitudes of the positive or negative peaks
%                   (depending on the value of the pos_neg parameter), cycle by cycle,
%                   of a voice segment
% i_P:               Positions of the positive or negative peaks (depending on the value
%                   of the parameter pos_neg), cycle by cycle, of a voice section;
%                   mPicosPos containing the peaks found in the analysis section
%                   (of the same duration as vSignal,).

% Check that the vector is of type column
if ~isvector( vSignal )
    error( 'vSignal is not a vector!' );
elseif size( vSignal, 2 ) ~= 1
    vSignal=vSignal';
end

P=[];
i_P=[];

Nsegm=length( vPitchPeriods );

% iFrame contains the durations in samples of the segments that has been used
% for the calculation of T (Each segment lasts 3 times the pitch period of the previous segment; the
% first segment lasts 34 ms)
iFrame=zeros( 1, Nsegm );
iFrame(1)=round( 0.034*iFs );

UltimoSegmentoSonoro=0;
for i=2:Nsegm
    if vPitchPeriods(i-1)~=0 
        iFrame(i)=3*vPitchPeriods(i-1);
        UltimoSegmentoSonoro=iFrame(i);
    elseif UltimoSegmentoSonoro~=0 && ~isnan( UltimoSegmentoSonoro )
        iFrame(i)=UltimoSegmentoSonoro;
    else
        iFrame(i)=iFrame(1);
    end
end

% For each segment the peaks are searched. There should be 3 peaks
% corresponding to the three pitch periods. Within each segment the point of
% start is the first zero crossing. From that point, and at an interval
% of T(n) samples, the first peak is searched. The following searches are made
% from the previous peak found, looking at distance T (n) samples, and in a
% 2RT(n) neighbourhood. T(n) is the pitch period for the segment in which we are located.
% If a segment is silent, it is skipped.

% k is an index to reference vectors P (peak amplitudes), and
% i_P (positions of these peaks on the complete signal s).
k=1;
R=0.2;

% Starting point of the segment analysed
IniSeg=iInicio;

% Reference sample
mueBusq=iInicio;

% comienzo=1 if we start the search from the first zero crossing (this happens if
% we are looking for the first peak), and is 0 if searching from the last peak found,
% to T(n) samples of distance in a 2RT(n) environment.
comienzo=1;

% For all segments
for n=1:Nsegm
    
    % If it is silent, the procedure is skiped
    if vPitchPeriods(n)==0 || isnan( vPitchPeriods(n) )
        IniSeg=IniSeg+iFrame(n);
        mueBusq=IniSeg;
        comienzo=1;
        
        % Otherwise, peaks are looked for
    else
        
        % While still in the segment
        while mueBusq+vPitchPeriods(n)<=IniSeg+iFrame(n)-1
            
            % comienzo = 1 if the first peak of the segment is looked for
            if comienzo==1
                
                % Starting point: first zero-crossing point
                if vSignal(mueBusq)==0
                    
                elseif vSignal(mueBusq)>0
                    while vSignal(mueBusq)>0 && (mueBusq<=IniSeg+iFrame(n)-1)
                        mueBusq=mueBusq+1;
                    end
                else
                    while vSignal(mueBusq)<0 && (mueBusq<=IniSeg+iFrame(n)-1)
                        mueBusq=mueBusq+1;
                    end
                end
                
                % The search of the first peak starts from the first zero-crossing
                if mueBusq+vPitchPeriods(n)<=IniSeg+iFrame(n)-1
                    if pos_neg==1
                        [P(k),pos]=max( vSignal(mueBusq:mueBusq+vPitchPeriods(n)) );
                    else
                        [P(k),pos]=min( vSignal(mueBusq:mueBusq+vPitchPeriods(n)) );
                    end
                    
                    i_P(k)=pos+(mueBusq-1);
                    % The peak serves as a reference for the next one
                    mueBusq=i_P(k);
                    k=k+1;
                    comienzo=0;
                end
                % If the first peak of the segment has already been found (start = 0)
            else
                % The peak found is the starting sample for the next search; from this point
                % the search occurs in T(n) samples, where n is the segment number where the peak is.
                margen=round(R*vPitchPeriods(n));
                mueIniBusq=mueBusq+vPitchPeriods(n)-margen;
                mueFinBusq=min(mueBusq+vPitchPeriods(n)+margen, IniSeg+iFrame(n)-1);
                try
                if pos_neg==1
                    [P(k),pos]=max( vSignal(mueIniBusq:mueFinBusq ));
                else
                    [P(k),pos]=min( vSignal(mueIniBusq:mueFinBusq ));
                end
                catch
                    clc
                end
                
                i_P(k)=pos+(mueIniBusq-1);
                % The current peak serves as a reference for the search of the
                % next peak
                mueBusq=i_P(k);
                k=k+1;
            end
            
        end
        
        % Next segment
        IniSeg=IniSeg+iFrame(n);
        
    end    
end

% Positive peaks
iNs=length( vSignal );
mPicosPos=zeros(1, iNs);
if ~isempty(P)
    mPicosPos(i_P)=P;
end

if nargout ==0
    eje_t=(iInicio:iFinal)/iFs;
    plot( eje_t, vSignal(iInicio:iFinal),'r', eje_t, mPicosPos(iInicio:iFinal),'g');
    if pos_neg==1
        title('positive peaks');
    else
        title('negative peaks');
    end
end