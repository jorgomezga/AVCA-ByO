function [P, i_P, picos]=PicosPos( vSignal, iFs, iInicio, iFinal, T )

% Gets the amplitudes P and the i_P positions of the positive peaks,
% cycle by cycle, of the input voice; it also returns a signal
% of the same duration of s, containing the peaks found
%
% Input parameters
%   T:           is the sequence of pitch period values (in samples)
%                calculated on windows of duration 3To without overlap,
%                where To is the pitch period in samples from the previous segment,
%                and 34ms the duration of the first window, as Boyanov does.
%                It is the result of applying the CalculatePitch function
%                using the Boyanov method
%
% Output parameters:
%   P:          amplitudes of positive peaks, cycle to cycle, of a speech
%               signal
%   i_P:        negative peak positions, cycle by cycle, of a speech signal
%   picos:      signal that contains the peaks found (of equal duration to vSignal)

if nargin<3, iInicio=1; end

P=[];
i_P=[];

Nsegm=length(T);

% The variable tam_seg contains the durations in samples of these segments that have been used
% for the calculation of T (Each segment lasts 3 times the pitch period of the previous segment; the
% first segment lasts 34 ms)
tam_seg=zeros(1,Nsegm);
tam_seg(1)=round(0.034*iFs);

UltimoSegmentoSonoro=0;
for i=2:Nsegm
    if T(i-1)~=0
        tam_seg(i)=3*T(i-1);
        UltimoSegmentoSonoro=tam_seg(i);
    elseif UltimoSegmentoSonoro~=0
        tam_seg(i)=UltimoSegmentoSonoro;
    else
        tam_seg(i)=tam_seg(1);
    end
end

% Positive peaks are searched for each segment. There should be 3 peaks
% corresponding to the three pitch periods. Within each segment the point of
% start of search is the first zero crossing. From that point, and at an interval
% of T(n) samples, the first positive peak is sought. The following searches are made
% from the previous positive peak found, looking at distance T(n) samples, and in a
% 2RT environment (n). T(n) is the pitch period for the segment in which we are located.
% If a segment is silent, it is skipped.

% k is an index to reference vectors P (amplitudes of positive peaks), and
% i_P (positions of these peaks on the complete signal s).
k=1;

R=0.2;

% Beginning
IniSeg=iInicio;

% Refence sample for the search
mueBusq=iInicio;

% start is 1 if we start the search from the first zero crossing (this occurs if
% we look for the first peak), and is 0 if we are searching from the last peak found,
% T(n) samples of distance in a 2RT(n) environment.
comienzo=1;

for n=1:Nsegm
    
    % Skip if the segment is silent
    if T(n)==0
        IniSeg=IniSeg+tam_seg(n);
        mueBusq=IniSeg;
        comienzo=1;
        
    else
        
        while mueBusq+T(n)<=IniSeg+tam_seg(n)-1
            
            % comienzo=1 for the first peak search
            if comienzo==1
                
                % Search for starting point: the first zero crossing
                if s(mueBusq)==0
                    
                elseif s(mueBusq)>0
                    while s(mueBusq)>0 && (mueBusq<=IniSeg+tam_seg(n)-1)
                        mueBusq=mueBusq+1;
                    end
                else
                    while s(mueBusq)<0 && (mueBusq<=IniSeg+tam_seg(n)-1)
                        mueBusq=mueBusq+1;
                    end
                end
                
                % The first peak from the first zero crossing to T (n) samples will be searched
                if mueBusq+T(n)<=IniSeg+tam_seg(n)-1
                    [P(k),pos]=max(s(mueBusq:mueBusq+T(n)));
                    i_P(k)=pos+(mueBusq-1);
                    % The first peak from the first zero crossing to T (n) samples will be searched
                    mueBusq=i_P(k);
                    k=k+1;
                    comienzo=0;
                end
                
                % If the first peak has been found (comienzo=0)
            else
                % The peak found is the starting sample for the next search; from that point
                % the search is done around T(n) samples, where n is the sample where the peak is.
                margen=round(R*T(n));
                mueIniBusq=mueBusq+T(n)-margen;
                mueFinBusq=min(mueBusq+T(n)+margen, IniSeg+tam_seg(n)-1);
                
                [P(k),pos]=max(s(mueIniBusq:mueFinBusq));
                i_P(k)=pos+(mueIniBusq-1);
                % The peak that is found serves as a reference for the next
                % search
                mueBusq=i_P(k);
                k=k+1;
            end
            
        end
        
        % We begin at the start of the next segment
        IniSeg=IniSeg+tam_seg(n);
        
    end
end

% Signal containing the positive peaks
Ns=length(s);
picos=zeros(1,Ns);
if ~isempty(P)
    picos(i_P)=P;
end

if nargout==0
    eje_t=(iInicio:length(vSignal))/iFs;
    plot(eje_t,vSignal(iInicio:end),'r',eje_t,picos(iInicio:end),'g');
    title('picos positivos');
end