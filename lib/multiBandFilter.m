function Y = multiBandFilter(signal,maxF,SR)

nChans  = size(signal,1);
nSamps  = size(signal,2);

[B,A] = singleHZBandPassFilters(maxF,SR);

Y = zeros(nChans,maxF,nSamps);
for ch = 1:nChans
    x = signal(ch,:);
    parfor ff = 1:maxF
        Y(ch,ff,:) = filtfilt(B(ff,:),A(ff,:),x);        
    end    
end