function [B,A] = singleHZBandPassFilters(maxFreq,fs)

Fr   = 1:maxFreq;
nF   = numel(Fr);

dP = 1;
dS = 40;
fOrd = 3;

B = zeros(nF,1+fOrd*2);
A = zeros(nF,1+fOrd*2);

for ff = 1:nF    
    [B(ff,:),A(ff,:)]=ellip(fOrd,dP,dS,(Fr(ff)+[-0.5 0.5])/fs*2);
end
