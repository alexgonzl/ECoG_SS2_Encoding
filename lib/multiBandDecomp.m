function [Amp,Ang] = multiBandDecomp(signal)

nChans  = size(signal,1);
Amp     = zeros(size(signal));
Ang     = zeros(size(signal));
for ch = 1:nChans
    Y = hilbert(signal(ch,:));
    Amp(ch,:) = abs(Y);
    Ang(ch,:) = angle(Y);
end