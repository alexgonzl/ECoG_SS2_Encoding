nB = out.nBands;
nC = out.nChans;
N = nB*nC;

bands = repmat((1:nBands)',[nChans,1]);
temp=repmat(1:nChans,[nBands,1]);
chans=temp(:);
    
rois = zeros(N,1);
for ii=1:3
    chanIDs = find(out.ROIid==ii);
    rois(ismember(chans,chanIDs))=ii;
end

bandroi = zeros(N,1);
cnt=1;
for ii=1:3
    for jj=1:6
        bandroi(bands==jj & rois== ii) = cnt;
        cnt = cnt+1;
    end
end
%%
x1=out.StudyRTs_PrePostActModel_Ts1(:,:,1); x1=x1(:);
x2=out.StudyRTs_PrePostActModel_Ts1(:,:,2); x2=x2(:);
y1=out.TestRTs_PrePostActModel_Ts1(:,:,1); y1=y1(:);
y2=out.TestRTs_PrePostActModel_Ts1(:,:,2); y2=y2(:);
%%
bandColors1 = brewermap(10,'Reds');
bandColors1(1:4,:) = [];
bandColors2 = brewermap(10,'Blues');
bandColors2(1:4,:) = [];
bandColors3 = brewermap(10,'Greens');
bandColors3(1:4,:) = [];
bandCols = [bandColors1;bandColors2;bandColors3];


s=scatterhist(x,y,'group',bandroi,'kernel','on','color',bandCols);

%%
%z=imagesc(corr(x1',x2'));
x1=out.StudyRTs_PrePostActModel_Ts1(:,:,1);
x2=out.StudyRTs_PrePostActModel_Ts1(:,:,2);
y1=out.TestRTs_PrePostActModel_Ts1(:,:,1); 
y2=out.TestRTs_PrePostActModel_Ts1(:,:,2); 
%%
r1=3;
r2=3;
[c,p]=corr(x1(:,data.ROIid==r1)',x2(:,data.ROIid==r2)');
imagesc(c.*(p<0.01),[-0.6 0.6])
%%
r1=2;
r2=3;
figure(1);
[c1,p]=corr(x1(:,data.ROIid==r1),x2(:,data.ROIid==r2));
imagesc(c1.*(p<0.01),[-0.6 0.6])

figure(2);
[c2,p]=corr(y1(:,data.ROIid==r1),y2(:,data.ROIid==r2));
imagesc(c2.*(p<0.01),[-0.6 0.6])

%%
chanIDByRegionSubj = cell(8,3);
chans = (1:134)';
for ss=1:8
    for jj = 1:3
        chanIDByRegionSubj{ss,jj} = chans(data.subjChans==ss & data.ROIid==jj);
    end
end

%% get ips/spl ips/ag and spl/ag chan pairs by subj


for ii = 1:3
    for jj=1:3
        
    end
end