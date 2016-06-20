
function exampleECoG_Sigs()
%% Load data
load('/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/raw_data/SRb/RawData/SRb-27/iEEGSRb-27_05.mat')
SR = 3051.76;
x = wave(:);
N = numel(x);
%% plot data segment.
tt = linspace(0,N/SR,N);
seed = 2;
rng(seed);
startPoint  = randi([1 round(N-SR)],1);
endPoint    = startPoint+ 2*round(SR);
epch        = startPoint:endPoint;
han = plotSignal(tt(epch),x(epch));

fPath = '/Users/alexandergonzalez/Google Drive/Research/Thesis/Defense/figures/';
print(han,'-dpdf', [fPath 'RawECoG_Epc'  num2str(seed)])
%% Spectrum
nfft = 2048;
window = hann(2048*2);

f = 0:0.5:500;
pxx = pwelch(x,window,numel(window)/2,f,SR,'power');

han = plotPSD(f,pxx);
pxx_yLims=ylim;
fPath = '/Users/alexandergonzalez/Google Drive/Research/Thesis/Defense/figures/';
print(han,'-dpdf', [fPath 'RawECoGSpec'])

%%
han = plotPSD(f,pxx);
p=patch([0.5 0.5 180 180],[pxx_yLims(1) pxx_yLims(2) pxx_yLims(2) pxx_yLims(1)],'r');
p.FaceAlpha=0.2;
p.EdgeColor='none';
fPath = '/Users/alexandergonzalez/Google Drive/Research/Thesis/Defense/figures/';
print(han,'-dpdf', [fPath 'RawECoGSpec2'])

%% Filter x:
lp = 180;
hp = 0.5;
notches = [60 120 180];
y = channelFilt(x',SR,lp,hp,[]);
for n = notches
    y = channelFilt(y,SR,[],[],n);
end
y = y(:);
%% Plot y
han = plotSignal(tt(epch),y(epch));
fPath = '/Users/alexandergonzalez/Google Drive/Research/Thesis/Defense/figures/';
print(han,'-dpdf', [fPath 'FiltECoG_Epc'  num2str(seed)])

han = plotSignal2(tt(epch),[x(epch) y(epch)],[0.5 0.5 0.5; 0.8 0 0]);
fPath = '/Users/alexandergonzalez/Google Drive/Research/Thesis/Defense/figures/';
print(han,'-dpdf', [fPath 'Filt2ECoG_Epc'  num2str(seed)])
%% spectrum of y
window = hann(2048*2);
f = 0:0.25:200;
pyy = pwelch(y,window,numel(window)/2,f,SR,'power');
%%
han = plotPSD(f,pyy);
ylim(pxx_yLims)
xlim([0 200]);
fPath = '/Users/alexandergonzalez/Google Drive/Research/Thesis/Defense/figures/';
print(han,'-dpdf', [fPath 'FiltECoGSpec'])

p=patch([80 80 180 180],[pxx_yLims(1) pxx_yLims(2) pxx_yLims(2) pxx_yLims(1)],'r');
p.FaceAlpha=0.2;
p.EdgeColor='none';
fPath = '/Users/alexandergonzalez/Google Drive/Research/Thesis/Defense/figures/';
print(han,'-dpdf', [fPath 'FiltECoGSpec2'])
%% color each band in the spectrum
han = plotPSD(f,pyy);
ylim(pxx_yLims)
xlim([0 200]);
bands         = [1 4;4 8;8 12;12 30; 30 80;80 180; 30 180;]; % in Hz
cols = brewermap(6,'Set3');
for jj = 1:6
    p=patch([bands(jj,1) bands(jj,1) bands(jj,2) bands(jj,2)],...
        [pxx_yLims(1) pxx_yLims(2) pxx_yLims(2) pxx_yLims(1)],'r');
    p.FaceColor = cols(jj,:);
    p.FaceAlpha =0.2;
    p.EdgeColor='none';
end
fPath = '/Users/alexandergonzalez/Google Drive/Research/Thesis/Defense/figures/';
print(han,'-dpdf', [fPath 'FiltECoGSpec3'])
%% color each band in the spectrum individually
bands         = [1 4;4 8;8 12;12 30; 30 80;80 180; 30 180;]; % in Hz
cols = brewermap(6,'Set3');
for jj = 1:6
    han = plotPSD(f,pyy);
ylim(pxx_yLims)
xlim([0 200]);

    p=patch([bands(jj,1) bands(jj,1) bands(jj,2) bands(jj,2)],...
        [pxx_yLims(1) pxx_yLims(2) pxx_yLims(2) pxx_yLims(1)],'r');
    p.FaceAlpha =0.2;
    p.EdgeColor='none';
    fPath = '/Users/alexandergonzalez/Google Drive/Research/Thesis/Defense/figures/';
    print(han,'-dpdf', [fPath 'FiltECoGSpec3_band' num2str(jj) ] )
end



% ylim([-40 25])
%% filter y to hgam
z = channelFilt(y',SR,180,80,[]); z = z(:);
han = plotSignal(tt(epch),z(epch));
fPath = '/Users/alexandergonzalez/Google Drive/Research/Thesis/Defense/figures/';
print(han,'-dpdf', [fPath 'FiltHGMECoG_Epc'  num2str(seed)])

%%
za = abs(hilbert(z));
han = plotSignal2(tt(epch),[z(epch) za(epch)],[0.5 0.5 0.5; 0.8 0 0]);
print(han,'-dpdf', [fPath 'HGAM_Hilb_ECoG_Epc'  num2str(seed)])

epch2 = epch(500:1200);
han = plotSignal2(tt(epch2),[z(epch2) za(epch2)],[0.5 0.5 0.5; 0.8 0 0]);
print(han,'-dpdf', [fPath 'HGAM_Hilb_Zommed_ECoG_Epc'  num2str(seed)])

zaf = channelFilt(za',SR,20,[],[]); zaf = zaf(:);
han = plotSignal2(tt(epch),[za(epch) zaf(epch)],[0.9 0.5 0.5; 0 0 0.5]);
print(han,'-dpdf', [fPath 'HGAM_Hilb2_ECoG_Epc'  num2str(seed)])

%
han = plotSignal2(tt(epch2),[za(epch2) zaf(epch2)],[0.9 0.5 0.5; 0 0 0.5]);
print(han,'-dpdf', [fPath 'HGAM_Hilb2_Zommed_ECoG_Epc'  num2str(seed)])

han = plotSignal(tt(epch2),20*log10(zaf(epch2)));
ylabel( '20log_{10}(|x|)')
print(han,'-dpdf', [fPath 'HGAM_Hilb_dB_Zommed_ECoG_Epc'  num2str(seed)])
%%
bands         = [1 4;4 8;8 12;12 30; 30 80;80 180; 30 180;]; % in Hz
bandNames    = {'delta','theta','alpha','beta','lgam','hgam'};

Z  = nan(6,N);
Z(6,:) = z;
Zaf =  nan(6,N);
Zaf(6,:) = zaf;
Zaf2 =  nan(6,N);
Zaf2(6,:) = zaf./mean(zaf);

for jj = 1:5
    Z(jj,:)  = channelFilt(y',SR,bands(jj,2),bands(jj,1),[]);
    Zaf(jj,:) = channelFilt(abs(hilbert(Z(jj,:))), SR, 20,[],[] );
    Zaf2(jj,:) = Zaf(jj,:)/mean(Zaf(jj,:));
end
%%
for jj=1:6
    han = plotSignal(tt(epch),20*log10(Zaf2(jj,epch)));
    ylabel('dB')
    %set(gca,'ytick',0:10:40)
    %ylim([-2 50])
    print(han,'-dpdf', [fPath bandNames{jj} 'Hilb2_ECoG_Epc'  num2str(seed)])
end
%%
freqs=1:6;
figure();
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[500 300],'paperposition',[0 0 500 300])
set(gcf,'position',[100,100,500,300]); % 100pt all around margin
a   = axes('position',[0.15 0.15 0.70 0.5]);
epch3 = epch(1:2:end);
Zafe = 20*log10((Zaf2(:,epch3)));
clev    = [-4:-1 1:4];
clev2 = quantile(Zafe(:),(clev+5)/10);
cmm       =     brewermap(numel(clev)-1,'*RdBu');
rcolmap = brewermap(1000,'*RdBu');
h=contourfcmap(tt(epch3),freqs, Zafe,clev2,cmm,'lo',rcolmap(1,:),'hi', rcolmap(end,:),'method','calccontour');
for i = 1:numel(h.l)
    h.l(i).Color='none';
end
hold on;
%plot([0 0],[1 6],'k','linewidth',4)
yticksLabels   = {'\delta','\theta','\alpha','\beta','l-\gamma','h-\gamma'};
set(gca,'fontsize',20,'ytick',freqs,'ytickLabel',yticksLabels)
print(gcf,'-dtiff','-r300', [fPath 'Spectrogram_ECoG_Epc'  num2str(seed)])

end
