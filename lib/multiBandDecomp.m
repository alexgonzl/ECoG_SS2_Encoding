function [Amp,Ang] = multiBandDecomp(signal)

N  = size(signal,1);
Amp     = zeros(size(signal));
Ang     = zeros(size(signal));
for ii = 1:N
    Y = hilbert(signal(ii,:));
    Amp(ii,:) = abs(Y);
    Ang(ii,:) = angle(Y);
end