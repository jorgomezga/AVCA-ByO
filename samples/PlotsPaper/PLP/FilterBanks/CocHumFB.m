  function H=CocHumFB(W)
% function H=CocHumFB(W);
%
% Filterbank responses according a cochlear model for the 
% human cochlea and for a samplig frequency of 8kHz. 
% (see CocHum.m)


%fs=8000;
%Nbank=35;

  load COCH
  H(1,:)=freqz(COCHB1,COCHA1,W);
   for k=1:Nbank-1,
    H(k+1,:)=H(k,:).*freqz(COCHB(k,:),COCHA(k,:),W);
   end

return


%--- phenomenological model (COCF)
%=================================
 
  load coc52
  h=freqz(Bm,Am,W);
  for k=1:Nbank,
    h=h.*freqz(Bt(k,:),At(k,:),W);
    hs=freqz(Bs(k,:),As(k,:),W);
    H(k,:)=h.*hs;
  end




