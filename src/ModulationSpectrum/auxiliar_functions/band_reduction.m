function [ vmMS_out ] = band_reduction(vmMS_in, iNBand  )

% [ vmMS_out ] = band_reduction(vmMS_in, iNBand  )
% Reduces the number of modulationbands in the input matrix vector
% (vmMs_in) to iNBand. The new output matrices are the result of the sum of
% squares of modulus of the original bands.
%
%   INPUTS:
%           -vmMs_in= Matrix vector containing Modulation spectra (one
%           modulation spectrum per matrix). For each matrix, rows index
%           refers to the modulation axis and the column index, to acoustic
%           axis.
%           - iNBand= Number of modulation bands of the new output
%           modulation spectra.
%   OUTPUTS:
%           - vmMs_out: Matrix vector containing the new modulation spectra
%           with iNBand moduation bands per matrix.
%

if nargin <2
    iNBand=32;
end

[iP,iQ,iR]=size(vmMS_in);

vmMS_out=zeros(iP,iNBand,iR);
nB=log2(iQ/2)/(ceil(iNBand/2));

vBandasp=power(2,(nB.*(0:ceil(iNBand/2))));

 %vBandas contains the frequency bands indices
vBandas=[ vBandasp, (iQ+1-flipdim(vBandasp(1:ceil(iNBand/2)),2))];
for j=1:iR
    vmMS_out(:,1,j)=vmMS_in(:,1,j);
    for i=1:iNBand
        % only for low frequency bands
        if (floor(vBandas(i+1))-ceil(vBandas(i)))<1 
            nA=0;
        else
            %sum of intermediate bands
            nA=sum(power(abs(vmMS_in(:,ceil(vBandas(i)+1):floor(vBandas(i+1)),j)),2),2); 
        end
        %lower residual band
        nB=(ceil(vBandas(i))-vBandas(i)).*power(abs(vmMS_in(:,ceil(vBandas(i)),j)),2);
        %upper residual band
        nC=(vBandas(i+1)-floor(vBandas(i+1))).*power(abs(vmMS_in(:,ceil(vBandas(i)),j)),2);
        
        vmMS_out(:,i,j)= sqrt(nA+nB+nC);
    end
end
        


end

