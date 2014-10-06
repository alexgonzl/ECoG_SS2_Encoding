function timeStamps = findEventMarkers(pDiodeSignal, matlabStamps)
% returns the event time stamps (in seconds) of the photo diode signal
% first input -> angl signal from photo diode
% second input -> maltab stimulus samples (VBLTimeStamp)

% defaults
PD_SR   = 24414.1;  % photo diode sampling rate
stimDist = 1.5;       % minimum distance allowed between stims in seconds

nEvents = numel(matlabStamps);

% correct signal for noise
tic
B     = fir1(300,0.01/PD_SR,'low');
pDiodeSignal = filtfilt(B,1,double(pDiodeSignal));
pDiodeSignal(pDiodeSignal>2.5)=5; pDiodeSignal(pDiodeSignal<=2.5)=0;
fprintf('Cleaning pdiode signal: time elapsed %g\n',toc)

% find indices of peaks that are at least Signal.stimDist*PD_SR apart (in
% samples)
[~,pDiodeStamps]=findpeaks(pDiodeSignal,'minpeakdistance',floor(stimDist*PD_SR),'minpeakheight',2.5);
pDiodeStamps = pDiodeStamps/PD_SR;

% match time derivatives
[c,lags]=crosscorr(diff(pDiodeStamps),diff(matlabStamps));
[~,id]=max(c);

timeStamps=lagmatrix(pDiodeStamps,lags(id)); 
timeStamps(isnan(timeStamps))=[]; 

if numel(timeStamps) < nEvents
    error('number of stamps do not match')
else
    timeStamps = timeStamps(1:nEvents);
end

assert(corr(timeStamps,matlabStamps')>0.999, 'bad match between diode and matlab samples')
fprintf('Match between diode and matlab is: %g%% \n',corr(timeStamps,matlabStamps')*100)

return