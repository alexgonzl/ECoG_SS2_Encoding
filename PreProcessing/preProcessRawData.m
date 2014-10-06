function data=preProcessRawData(subj,expt)
% function filters the raw data traces and concantenates data blocks for
% later processing
%
%   inputs:
%       subuject name and experiments
%
%   data dependencies:
%       photo-diode event markers in the same sample space as the data
%       Raw data for each block
%       
%   file dependencies: 
%       subjChanInfo.m
%       subjExptInfo.m
%       channelFilt.m
%       ts_adaptivezscore.m
% 
%   Author: A.Gonzl. 
%% intial settings

data                = [];
data.subjNum        = subj;

temp                = subjChanInfo(data.subjNum);
data.chanInfo       = temp;
data.nChans         = numel(data.chanInfo.allChans);

temp                = SS2e_subjInfo(data,subj);
data.subjName       = temp.subjName;
data.blocklist      = temp.blocklist;   blocklist = data.blocklist;
data.nBlocks        = numel(data.blocklist);

data.paths.dataPath         = ['/Volumes/ECoG_SS2/SS2/data/' data.subjName '/' ];
data.paths.behavDataPath    = [data.paths.dataPath 'BehavData/'];
data.paths.RawDataPath      = [data.paths.dataPath 'RawData/'];
data.paths.ResultsPath      = ['/Volumes/ECoG_SS2/SS2/SS2e/Results/' data.subjName '/' ];
data.paths.preProRawPath    = [data.paths.ResultsPath 'preProcessed/'];         

load([data.paths.ResultsPath 'BehavResults/behavResults.mat'])
data.behavorInfo    = S;
data.eventTimeStamps = S.eventTimeStamps;
clear S;

data.standarized    = false;
if data.standarized, zscore_str = 'zscored';else zscore_str='';end
%temp = subjExptInfo(subj,expt,'allChCAR');

%% Parameters

data.SRorig     = 3051.76;
data.lowpass    = 180;  lp = data.lowpass; % low pass filter
data.hipass     = 1;    hp = data.hipass;% high pass filter
data.notch      = [60 120]; notches = data.notch; % notch; notches =
data.comp       = 7; comp=data.comp;% compression factor
data.SR         = data.SRorig/data.comp; SR = data.SR;

%% Decompose channels and save
rawDataPath             = data.paths.RawDataPath;
data.trialOnsets        = [];
data.blockOffset(1)     = 0;
refChan = data.chanInfo.refChannel; data.refChan = refChan;

data.signal =[];
for b = 1: data.nBlocks
    
    display(['Prossesing Block ' data.blocklist{b} ])
    
    % event time points
    stamps          = data.eventTimeStamps{b}; % in seconds.
    data.nevents(b) = numel(stamps);
    
    % pre-allocation
    rawdatafile = [rawDataPath data.blocklist{b} '/iEEG' data.blocklist{b} '_01.mat'];
    
    x           = load(rawdatafile);
    nBlockSamps = ceil(length(x.wave)/data.comp); clear x;
    
    data.blockOffset(b+1) = data.blockOffset(b)+nBlockSamps;
    data.trialOnsets = [data.trialOnsets ; ceil(stamps*data.SR)+data.blockOffset(b)];
    
    signal = zeros(data.nChans,nBlockSamps);
        
    for ch = 1:data.nChans    
        if (ch~=refChan)            
            channel = num2str(ch);
            display(['Prossesing Channel ' channel])
            zerosstr = num2str(zeros(1,2-numel(channel)));
            rawdatafile = [rawDataPath blocklist{b} ...
                '/iEEG' blocklist{b} '_' zerosstr channel '.mat'];
            
            % load data
            x = load(rawdatafile);
            x = -double(x.wave);
            % detrend, downsample, bandpass & notch
            x = detrend(x,'linear');
            x = decimate(x,comp);
            x = channelFilt(x,SR,lp,hp,[]);
            for n = notches
                x = channelFilt(x,SR,[],[],n);
            end            
            if data.standarized
                x = ts_adaptivezscore(x);
            end
            signal(ch,:) = x;
        end
    end
    data.signal = cat(2,data.signal,signal);
end

if ~exist(data.paths.preProRawPath,'dir'), mkdir(data.paths.preProRawPath),end;
save([data.paths.preProRawPath '/data' zscore_str '.mat'] , 'data')

