function [result1] = measuremnet_types(xr1, xr2, model)
%xr = zscore(xr);%z-score normalize each column of matrix xr
HT1 = hilbert(xr1);
HT2 = hilbert(xr2);

%ph = atan( imag(HT)./real(HT) );
ph2 = angle(HT2);
ph1 = angle(HT1);


if model=="mod(π)_PLV"
    phi = mod(unwrap(ph1),pi) - mod(unwrap(ph2),pi);
    E = exp(1i*phi);
    result1 = abs(sum(E,'all'))/size(xr1,2);
elseif model=="PLV"
    phi = ph1 - ph2;
    E = exp(1i*phi);
    result1 = abs(sum(E,'all'))/size(xr1,2);
elseif model=="mod(π)_Phase_Diff"
    phi = mod(ph1,pi) - mod(ph2,pi);
    E = (1-(abs(phi/pi)));
    result1 = sum(E,'all')/size(xr1,1);
elseif model=="Std_Diff"
    xr1 =(xr1 - min(xr1)) / ( max(xr1) - min(xr1) );
    xr2 =(xr2 - min(xr2)) / ( max(xr2) - min(xr2) );
    result1 = 1-std(xr1-xr2);
elseif model=="LoS"
    result1 = sum(gaussmf(xr1 - xr2,[0.1 0]),'all')/size(xr1,2);
end
