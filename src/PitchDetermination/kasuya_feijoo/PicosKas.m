function [P, i_P, picos]=PicosKas( vSignal, iFs, iInicio, iFinal, T, iPosNeg )

%       Computes the amplitudes P and the i_P positions of the positive or negative peaks
%       (depending on the iPosNeg parameter), cycle by cycle, of a speech
%       section; it also returns the signal peaks of the same duration as vSignal, i.e. the peaks
%       found in the analysis.
%
% Input parameters:
%   vSignal: is a column vector containing the complete speech signal
%   iFs:     is the sampling frequency
%   iIni:    is the first sample of the section to be analyzed
%   iFinal:  is the last sample of the section to be analyzed
%   T:       is the sequence of pitch period values (in samples)
%            calculated on 40ms windows, with displacements of 20ms,
%            as Kasuya and Feijoo do. 
%       Is the result of applying the PitchKas function.
%   iPosNeg: indicates the peaks to be searched. 
%           If iPosNeg = 1 -> positive.
%           If iPosNeg = -1: negative
%
% Output parameters:
% P:
% i_P:
% peaks:

% Check that the vector is of type column
if ~isvector( vSignal ) 
    error( 'vSignal is not a vector!' );
elseif size( vSignal, 2 ) ~= 1 
    vSignal=vSignal'; 
end

Nsegm=length(T);

% The search for peaks is done throughout the analysis section.
% The starting point is the first zero crossing. From that point, and at an interval
% of T(1) samples the first peak is searched. The next searches are made starting from the
% previously found peak, within a distance of T(n) samples, and in a 2RT(n) environment, R = 0.2.
% T(n) is the pitch period for the segment in which we are located;
% as there is overlap of windows we consider that: segment(1)= zona_sin_solape_ventana(1),
% segment(2)=zona_solapada_ventanas (1) and (2), etc. If a segment is silent the search is skipped.
P=[];
i_P=[];
sin_solape=0.02*iFs;

% k is an index to reference vectors P (peak amplitudes), and
% i_P (positions of these peaks on the complete vSignal signal).
k = 1;

% n is used to indicate the segment we are in
n = 1;

R=0.2;

% Reference sample
mueBusq=iInicio;

% start is 1 if we start the search from the first zero crossing (this happens if
% we are looking for the first peak, or if we just skipped a silent segment), and it is 0 if
% the search starts from the last peak found, within a distance of T(n) samples in a 2RT (n) environment.
comienzo=1;
terminado=0;

% If silent segments, skip
while (T(n)==0) && (n<Nsegm)
   n=n+1;
   mueBusq=mueBusq+sin_solape;
end

if (n==Nsegm) && T(n)==0
   terminado=1;
end


while (mueBusq+T(n)<=iFinal) && (terminado==0)
   
   if comienzo==1
       % Search of the starting point
       if vSignal(mueBusq)==0
           
      elseif vSignal(mueBusq)>0
         while vSignal(mueBusq)>0 && (mueBusq<=iFinal)
            mueBusq=mueBusq+1;
         end
      else
         while vSignal(mueBusq)<0 && (mueBusq<=iFinal)
            mueBusq=mueBusq+1;
         end
      end
   
      % Search of the first peak (or the first one after a silent segment)
      if mueBusq+T(n)<=iFinal
         if iPosNeg==1
            [P(k),pos]=max( vSignal(mueBusq:mueBusq+T(n)) );
         else
            [P(k),pos]=min( vSignal(mueBusq:mueBusq+T(n)) );           
         end
         
         i_P(k)=pos+(mueBusq-1);
         % The found peak, serves as a reference for the next
         mueBusq=i_P(k);
         k=k+1;
         comienzo=0;
      end
      
   else 
      % The peak found is the starting (iInicio) sample for the next
      % search; from there, it is looked for at T(n) samples, where n is 
      % the sample number where the peak is.
      margen=round(R*T(n));
      mueIniBusq=mueBusq+T(n)-margen;
      mueFinBusq=min(mueBusq+T(n)+margen, iFinal);
      
      if iPosNeg==1
         [P(k),pos]=max(vSignal(mueIniBusq:mueFinBusq));
      else
         [P(k),pos]=min(vSignal(mueIniBusq:mueFinBusq));
      end
      
      i_P(k)=pos+(mueIniBusq-1);      
      % The found peak, serves as a reference for the next
      mueBusq=i_P(k);
      k=k+1;

      n=ceil( (mueBusq-iInicio)/sin_solape );
      if n>Nsegm
         n=Nsegm;
      end
         
      % If the segment is silent, then skip
      while (T(n)==0) && (n<Nsegm)    
         comienzo=1;
         n=n+1;
         mueBusq=mueBusq+sin_solape;
      end
         
      if (n==Nsegm) && T(n)==0
         terminado=1;
      end      
   end
   
end

% Signal containing the peaks
Ns=length(vSignal);
picos=zeros(1,Ns);
if ~isempty(P)
   picos(i_P)=P;
end

if nargout == 0 
      eje_t=(iInicio:iFinal)/iFs;
      plot( eje_t, vSignal(iInicio:iFinal), 'r', eje_t, picos(iInicio:iFinal), 'g');
      if iPosNeg==1
          title('positive peaks');
      else
          title('negative peaks'); 
      end
end