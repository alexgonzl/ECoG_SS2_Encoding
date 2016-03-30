function PCAres=MMPCA(Y,R)

% performs PCA based on the multivariate MM estimate of shape, with
% fast and robust bootstrap
%
% [calls multiMM.m and multiMMRobBoot.m]
%
% INPUT :
%   Y : (n x q) data
%   R : number of bootstrap samples
% OUTPUT :
%   PCAres.MMest : (structure) result of multiMM
%   PCAres.MMboot : (structure) result of multiMMRobBoot
%   PCAres.vecMMest : (2*q+2*q*q x 1) MM(and S) estimates in vec-form
%   PCAres.eigG : (q x 1) eigenvalues of MM shape in ascending (!) order
%   PCAres.eigvec : (q*q x 1) eigenvectors of MM-shape (in vec-form)
%   PCAres.pvar : (q-1 x 1) percentages of variance for MM eigenvalues
%   PCAres.variance_MM : (2*q+2*q*q x 1) bootstrap variance for
%                                                       MM-estimates
%   PCAres.variance_eigG : (q x 1) bootstrap variance for MM eigenvalues
%   PCAres.variance_eigvec : (q*q x 1) bootstrap variance for MM eigenvectors
%   PCAres.variance_pvar : (q-1 x 1) bootstrap variance for percentage of
%                                               variance for MM-eigenvalues
%   PCAres.avgangle : (q x 1) average angles between bootstrap eigenvectors
%                                           and original MM eigenvectors
%   PCAres.CI95_MM : (2*q+2*q*q x 2) 95% BCa intervals for MM-estimates
%                                                           ([lower upper])
%   PCAres.CI95_eigG : (q x 2) 95% BCa intervals for MM eigenvalues
%   PCAres.CI95_eigvec : (q*q x 2) 95% BCa intervals for MM eigenvectors
%   PCAres.CI95_pvar : (q-1 x 2) 95% BCa intervals for percentage of
%                                               variance for MM-eigenvalues
%   PCAres.CI95one_pvar : (q-1 x 1) 95% one-sided BCa intervals for 
%               percentage of variance for MM-eigenvalues ([-infty upper]) 
%   PCAres.failedsamples : number of discarded bootstrap samples due to
%                                   non-positive definiteness of shape
%                   (results from 'discarded' bootstrap samples are omitted
%            such that actual number of recalculations can be lower than R) 
%   REMARK: whenever matrices contain eigenvalue-related values, the ordering is "ASCENDING"
%                       (e.g. first column of eigenvectors contains vector corresponding to smallest eigenvalue)                
%                       (or first row in case of the matrix of bootstrapped values)

[n q] = size(Y);
bdp = .5;

% compute the MM-estimates
MMests=multiMM(Y);

MMloc = MMests.loc;
MMSigma = MMests.covariance;
MMGamma = MMests.shape;
SSigma = MMests.Sshape*(MMests.auxscale)^2;
Sloc = MMests.Sloc;

% compute eigenvalues/vectors
[VGamma DGamma] = eig(MMGamma);
[VSigma DSigma] = eig(MMSigma);
[MMeigvalsGamma IXGamma] = sort(diag(DGamma));
[MMeigvalsSigma IXSigma] = sort(diag(DSigma));
for k=1:q
    [dummy IX]=sort(abs(VSigma(:,IXSigma(k))));
    VSigma(:,IXSigma(k))=sign(VSigma(IX(q),IXSigma(k)))*VSigma(:,IXSigma(k));
end    
MMeigvecs=vecop(VSigma(:,IXSigma));

MMvarperc=zeros(1,q-1);
for k=1:(q-1)
    MMvarperc(k)=sum(MMeigvalsSigma((q-k+1):q))/sum(MMeigvalsSigma);
end

% perform fast and robust bootstrap on MM-result
bootres = multiMMRobBoot(Y, ones(n,1), Sloc, SSigma, MMloc, MMSigma, R, bdp);

avgangle = mean(bootres.angles,2)';

MMvariances = var(bootres.centered,0,2)';
eigGvariances = var(bootres.eigGcentered,0,2)';
pvarvariances = var(bootres.pvarcentered,0,2)';
eigvecvariances = var(bootres.eigveccentered,0,2)';

% sort bootstrap recalculations for constructing intervals
sortedMMest = sort(bootres.centered,2);
sortedeigG = sort(bootres.eigGcentered,2);
sortedpvar = sort(bootres.pvarcentered,2);
sortedeigvec = sort(bootres.eigveccentered,2);

vecMMest = bootres.MMest;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute BCa confidence limits

% empirical inlfuences for computing a in BCa intervals, based on IF(MM)
EinfMM = einfsMM(Y, MMloc, MMSigma);
inflE = [EinfMM.loc EinfMM.shape EinfMM.covS EinfMM.locS];
eigGinflE = EinfMM.eigs;
eigvecinflE = EinfMM.eigvec;
pvarinflE = EinfMM.varperc;

% set R equal to the actual number of OK bootstrap samples
failedsamples = R - size(bootres.centered,2);
R = R - failedsamples;
MMest_95CIbca = zeros(2*q+2*q*q,2);
for i=1:(2*q+2*q*q)
    nofless=length(sortedMMest(i,sortedMMest(i,:)<=0));
    w=norminv(nofless/(R+1));
    a=1/6*sum(inflE(:,i).^3)/(sum(inflE(:,i).^2)^(3/2));
    alphatildelow95=normcdf(w+(w-1.96)/(1-a*(w-1.96)));
    alphatildehigh95=normcdf(w+(w+1.96)/(1-a*(w+1.96)));
    indexlow95=max((R+1)*alphatildelow95,1);
    indexlow95=min(indexlow95,R);
    indexhigh95=min((R+1)*alphatildehigh95,R);
    indexhigh95=max(indexhigh95,1);
    MMest_95CIbca(i,1) = sortedMMest(i,round(indexlow95)) + vecMMest(i);
    MMest_95CIbca(i,2) = sortedMMest(i,round(indexhigh95)) + vecMMest(i);
end
eigG_95CIbca = zeros(q,2);
for i=1:q
    nofless=length(sortedeigG(i,sortedeigG(i,:)<=0));
    w=norminv(nofless/(R+1));
    a=1/6*sum(eigGinflE(:,i).^3)/(sum(eigGinflE(:,i).^2)^(3/2));
    alphatildelow95=normcdf(w+(w-1.96)/(1-a*(w-1.96)));
    alphatildehigh95=normcdf(w+(w+1.96)/(1-a*(w+1.96)));
    indexlow95=max((R+1)*alphatildelow95,1);
    indexlow95=min(indexlow95,R);
    indexhigh95=min((R+1)*alphatildehigh95,R);
    indexhigh95=max(indexhigh95,1);
    eigG_95CIbca(i,1) = sortedeigG(i,round(indexlow95)) + MMeigvalsGamma(i);
    eigG_95CIbca(i,2) = sortedeigG(i,round(indexhigh95)) + MMeigvalsGamma(i);
end
eigvec_95CIbca = zeros(q*q,2);
for i=1:(q*q)
    nofless=length(sortedeigvec(i,sortedeigvec(i,:)<=0));
    w=norminv(nofless/(R+1));
    a=1/6*sum(eigvecinflE(:,i).^3)/(sum(eigvecinflE(:,i).^2)^(3/2));
    alphatildelow95=normcdf(w+(w-1.96)/(1-a*(w-1.96)));
    alphatildehigh95=normcdf(w+(w+1.96)/(1-a*(w+1.96)));
    indexlow95=max((R+1)*alphatildelow95,1);
    indexlow95=min(indexlow95,R);
    indexhigh95=min((R+1)*alphatildehigh95,R);
    indexhigh95=max(indexhigh95,1);
    eigvec_95CIbca(i,1) = sortedeigvec(i,round(indexlow95)) + MMeigvecs(i);
    eigvec_95CIbca(i,2) = sortedeigvec(i,round(indexhigh95)) + MMeigvecs(i);
end
pvar_95CIbca = zeros(q-1,2);
pvar_95CIbca_one = zeros(q-1,1);
for i=1:(q-1)
    nofless=length(sortedpvar(i,sortedpvar(i,:)<=0));
    w=norminv(nofless/(R+1));
    a=1/6*sum(pvarinflE(:,i).^3)/(sum(pvarinflE(:,i).^2)^(3/2));
    alphatildelow95=normcdf(w+(w-1.96)/(1-a*(w-1.96)));
    alphatildehigh95=normcdf(w+(w+1.96)/(1-a*(w+1.96)));
    indexlow95=max((R+1)*alphatildelow95,1);
    indexlow95=min(indexlow95,R);
    indexhigh95=min((R+1)*alphatildehigh95,R);
    indexhigh95=max(indexhigh95,1);
    pvar_95CIbca(i,1) = sortedpvar(i,round(indexlow95)) + MMvarperc(i);
    pvar_95CIbca(i,2) = sortedpvar(i,round(indexhigh95)) + MMvarperc(i);

    % one-sided interval
    alphatildehigh95=normcdf(w+(w+1.645)/(1-a*(w+1.645)));
    indexhigh95=min((R+1)*alphatildehigh95,R);
    indexhigh95=max(indexhigh95,1);
    pvar_95CIbca_one(i) = sortedpvar(i,round(indexhigh95)) + MMvarperc(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PCAres.MMest = MMests;
PCAres.MMboot = bootres;

PCAres.vecMMest = vecMMest;
PCAres.eigG = MMeigvalsGamma';
PCAres.eigvec = MMeigvecs;
PCAres.pvar = MMvarperc';

PCAres.variance_MM = MMvariances;
PCAres.variance_eigG = eigGvariances;
PCAres.variance_eigvec = eigvecvariances;
PCAres.variance_pvar = pvarvariances;
PCAres.avgangle = avgangle;

PCAres.CI95_MM = MMest_95CIbca;
PCAres.CI95_eigG = eigG_95CIbca;
PCAres.CI95_eigvec = eigvec_95CIbca;
PCAres.CI95_pvar = pvar_95CIbca;

PCAres.CI95one_pvar = pvar_95CIbca_one;

PCAres.failedsamples = failedsamples;

% --------------------------------------------------------------------

function einfs = einfsMM(Y, MMloc, MMcov)

% empirical influences for MM-estimates + PCA estimates

[n p] = size(Y);

c = biweight95eff(p);
[c0 b0] = multipar(p,0.5);

MMshape = det(MMcov)^(-1/p) * MMcov;

% compute eigenvalues/vectors
[VGamma DGamma] = eig(MMshape);
[eigs IXGamma] = sort(diag(DGamma));
for k = 1:p
    [dummy IX] = sort(abs(VGamma(:,IXGamma(k))));
    VGamma(:,IXGamma(k)) = sign(VGamma(IX(p),IXGamma(k)))*VGamma(:,IXGamma(k));
end    
eigvecs = VGamma(:,IXGamma);

varperc = zeros(1,p-1);
for k = 1:(p-1)
    varperc(k) = sum(eigs((p-k+1):p))/sum(eigs);
end

divec = sqrt(mahsq(Y, MMloc, MMcov, n));

psidervec = rhobiweightder2(divec, c);
psidervecS = rhobiweightder2(divec, c0);
uvec = rhobiweightder1(divec, c) ./ divec;
uvecS = rhobiweightder1(divec, c0) ./ divec;
vvec = rhobiweightder1(divec,c) .* divec;
vvecS = rhobiweightder1(divec,c0) .* divec;
rhovecS = rhobiweight(divec,c0);

betaMM = (1-1/p) * mean(uvec) + 1/p * mean(psidervec);
einfs.loc = 1 / betaMM * (uvec' * ones(1,p)) .* (Y - ones(n, 1) * MMloc);

betaS = (1-1/p) * mean(uvecS) + 1/p * mean(psidervecS);
einfs.locS = 1 / betaS * (uvecS' * ones(1,p)) .* (Y - ones(n, 1) * MMloc);

gamma3 = mean(vvecS);
gamma1MM = mean(psidervec.*(divec.^2) + (p+1)*vvec) / (p+2);
gamma1S = mean(psidervecS.*(divec.^2) + (p+1)*vvecS) / (p+2);

einfscov = zeros(n,p*p);
einfscovS = zeros(n,p*p);
einfsshape = zeros(n,p*p);
einfseigvec = zeros(n,p*p);
einfseigscov = zeros(n,p);
einfseigs = zeros(n,p);
einfsvarperc = zeros(n,p-1);
for i = 1:n
    IFcov = 2/gamma3*(rhovecS(i) - b0) * MMcov + 1/gamma1MM*p*vvec(i)*((Y(i,:)-MMloc)'*(Y(i,:)-MMloc)/divec(i)^2-1/p*MMcov);
    einfscov(i,:) = vecop( IFcov )';
    IFcovS = 2/gamma3*(rhovecS(i) - b0) * MMcov + 1/gamma1S*p*vvecS(i)*((Y(i,:)-MMloc)'*(Y(i,:)-MMloc)/divec(i)^2-1/p*MMcov);
    einfscovS(i,:) = vecop( IFcovS )';
    IFshape = 1/gamma1MM*p*vvec(i)*det(MMcov)^(-1/p)*((Y(i,:)-MMloc)'*(Y(i,:)-MMloc)/divec(i)^2-1/p*MMcov);
    einfsshape(i,:) = vecop( IFshape )';
    IFeigvecs = zeros(p,p);
    for k=1:p
        IFeigveck = zeros(p,1);
        for j=1:p
            if j~=k
                IFeigveck = IFeigveck + 1/(eigs(k)-eigs(j))*(eigvecs(:,j)'*IFshape*eigvecs(:,k))*eigvecs(:,j);
            end
        end
        IFeigvecs(:,k) = IFeigveck;
        einfseigs(i,k) = eigvecs(:,k)' * IFshape * eigvecs(:,k);
        einfseigscov(i,k) = eigvecs(:,k)' * IFcov * eigvecs(:,k);
    end
    einfseigvec(i,:) = vecop(IFeigvecs)';
    for k=1:(p-1)
        einfsvarperc(i,k) = 1/sum(eigs) * ((1-varperc(k))*sum(einfseigs(i,(p-k+1):p)) - varperc(k)*sum(einfseigs(i,1:(p-k))));
    end
end
einfs.shape = einfsshape;
einfs.cov = einfscov;
einfs.covS = einfscovS;
einfs.eigvec = einfseigvec;
einfs.eigs = einfseigs;
einfs.eigscov = einfseigscov;
einfs.varperc = einfsvarperc;

% --------------------------------------------------------------------

function rho=rhobiweight(x,c)

% Computes Tukey's biweight rho function met constante c voor alle waarden
% in de vector x.

hulp=x.^2/2-x.^4/(2*c^2)+x.^6/(6*c^4);
rho=hulp.*(abs(x)<c)+c^2/6.*(abs(x)>=c);

% --------------------------------------------------------------------

function psi=rhobiweightder1(x,c)

% Computes Tukey's biweight psi function met constante c voor alle waarden
% in de vector x.

hulp=x-2.*x.^3/(c^2)+x.^5/(c^4);
psi=hulp.*(abs(x)<c);

% --------------------------------------------------------------------

function psi=rhobiweightder2(x,c)

% Computes Tukey's biweight psi function met constante c voor alle waarden
% in de vector x.

hulp=1-6.*x.^2/(c^2)+5.*x.^4/(c^4);
psi=hulp.*(abs(x)<c);

% --------------------------------------------------------------------
 
function mah=mahsq(dat,meanvct,covmat,n)

% Computes the mahalanobis distances.

invcov=inv(covmat);

hlp=dat-ones(n,1)*meanvct;
mah=sum(hlp*invcov.*hlp,2)';

% --------------------------------------------------------------------

function vec=vecop(mat)
% performs vec-operation (stacks colums of a matrix into column-vector)

[nr nc]=size(mat);
vecmat=zeros(nr*nc,1);
for col=1:nc
    startindex=(col-1)*nr+1;
    vecmat(startindex:(startindex+nr-1))=mat(:,col);
end

vec=vecmat;

% --------------------------------------------------------------------

function c=biweight95eff(p)

x=[1    4.6851
   2    5.1230
   3    5.4903
   4    5.8103
   5    6.0963
   6    6.3562
   7    6.5956
   8    6.8182
   9    7.0268
   10   7.2235];
c=x(p,2);

% --------------------------------------------------------------------

function [c,b]=multipar(p,bdp)

x1=[1.0000    4.0963    0.4195
    2.0000    5.9815    0.8944
    3.0000    7.3996    1.3689
    4.0000    8.5863    1.8431
    5.0000    9.6277    2.3173
    6.0000   10.5669    2.7915
    7.0000   11.4291    3.2656
    8.0000   12.2307    3.7398
    9.0000   12.9829    4.2139
   10.0000   13.6938    4.6880
   11.0000   14.3696    5.1621
   12.0000   15.0150    5.6362
   13.0000   15.6337    6.1103
   14.0000   16.2289    6.5844
   15.0000   16.8030    7.0585
   16.0000   17.3582    7.5327
   17.0000   17.8961    8.0068
   18.0000   18.4183    8.4809
   19.0000   18.9261    8.9550
   20.0000   19.4207    9.4291
   21.0000   19.9029    9.9032
   22.0000   20.3738   10.3773
   23.0000   20.8340   10.8514
   24.0000   21.2842   11.3255
   25.0000   21.7252   11.7996
   26.0000   22.1573   12.2737
   27.0000   22.5812   12.7478
   28.0000   22.9973   13.2219
   29.0000   23.4060   13.6960
   30.0000   23.8076   14.1701
   31.0000   24.2026   14.6442
   32.0000   24.5913   15.1183
   33.0000   24.9739   15.5924
   34.0000   25.3507   16.0665
   35.0000   25.7220   16.5406
   36.0000   26.0881   17.0147
   37.0000   26.4490   17.4888
   38.0000   26.8051   17.9629
   39.0000   27.1566   18.4370
   40.0000   27.5035   18.9111
   41.0000   27.8461   19.3852
   42.0000   28.1846   19.8593
   43.0000   28.5190   20.3334
   44.0000   28.8496   20.8075
   45.0000   29.1764   21.2816
   46.0000   29.4996   21.7557
   47.0000   29.8193   22.2298
   48.0000   30.1356   22.7039
   49.0000   30.4486   23.1780
   50.0000   30.7585   23.6521
   51.0000   31.0652   24.1262];


x2=[1.0000    2.9370    0.3594
    2.0000    4.4274    0.8168
    3.0000    5.5281    1.2733
    4.0000    6.4426    1.7295
    5.0000    7.2423    2.1854
    6.0000    7.9619    2.6413
    7.0000    8.6215    3.0971
    8.0000    9.2342    3.5529
    9.0000    9.8086    4.0087
   10.0000   10.3511    4.4644
   11.0000   10.8666    4.9201
   12.0000   11.3587    5.3759
   13.0000   11.8304    5.8316
   14.0000   12.2839    6.2873
   15.0000   12.7213    6.7430
   16.0000   13.1441    7.1987
   17.0000   13.5538    7.6544
   18.0000   13.9514    8.1101
   19.0000   14.3380    8.5658
   20.0000   14.7144    9.0214
   21.0000   15.0815    9.4771
   22.0000   15.4398    9.9328
   23.0000   15.7900   10.3885
   24.0000   16.1326   10.8442
   25.0000   16.4681   11.2999
   26.0000   16.7968   11.7556
   27.0000   17.1193   12.2112
   28.0000   17.4358   12.6669
   29.0000   17.7466   13.1226
   30.0000   18.0521   13.5783
   31.0000   18.3525   14.0340
   32.0000   18.6481   14.4896
   33.0000   18.9391   14.9453
   34.0000   19.2256   15.4010
   35.0000   19.5080   15.8567
   36.0000   19.7863   16.3124
   37.0000   20.0607   16.7680
   38.0000   20.3315   17.2237
   39.0000   20.5987   17.6794
   40.0000   20.8624   18.1351
   41.0000   21.1229   18.5908
   42.0000   21.3802   19.0464
   43.0000   21.6345   19.5021
   44.0000   21.8858   19.9578
   45.0000   22.1342   20.4135
   46.0000   22.3799   20.8691
   47.0000   22.6229   21.3248
   48.0000   22.8633   21.7805
   49.0000   23.1013   22.2362
   50.0000   23.3368   22.6918
   51.0000   23.5699   23.1475];

x3=[1.0000    1.5476    0.1996
    2.0000    2.6608    0.5900
    3.0000    3.4529    0.9935
    4.0000    4.0966    1.3985
    5.0000    4.6520    1.8034
    6.0000    5.1477    2.2082
    7.0000    5.5995    2.6128
    8.0000    6.0173    3.0173
    9.0000    6.4078    3.4217
   10.0000    6.7758    3.8260
   11.0000    7.1248    4.2302
   12.0000    7.4574    4.6344
   13.0000    7.7758    5.0386
   14.0000    8.0816    5.4427
   15.0000    8.3763    5.8468
   16.0000    8.6609    6.2509
   17.0000    8.9364    6.6550
   18.0000    9.2037    7.0590
   19.0000    9.4634    7.4631
   20.0000    9.7162    7.8671
   21.0000    9.9626    8.2711
   22.0000   10.2030    8.6751
   23.0000   10.4379    9.0792
   24.0000   10.6676    9.4832
   25.0000   10.8925    9.8872
   26.0000   11.1128   10.2912
   27.0000   11.3288   10.6952
   28.0000   11.5408   11.0991
   29.0000   11.7489   11.5031
   30.0000   11.9535   11.9071
   31.0000   12.1546   12.3111
   32.0000   12.3524   12.7151
   33.0000   12.5471   13.1191
   34.0000   12.7388   13.5230
   35.0000   12.9276   13.9270
   36.0000   13.1138   14.3310
   37.0000   13.2973   14.7350
   38.0000   13.4784   15.1389
   39.0000   13.6570   15.5429
   40.0000   13.8334   15.9469
   41.0000   14.0075   16.3508
   42.0000   14.1795   16.7548
   43.0000   14.3494   17.1587
   44.0000   14.5173   17.5627
   45.0000   14.6833   17.9667
   46.0000   14.8475   18.3706
   47.0000   15.0098   18.7746
   48.0000   15.1705   19.1785
   49.0000   15.3294   19.5825
   50.0000   15.4867   19.9865
   51.0000   15.6424   20.3904];

if bdp==0.5
    c=x3(p,2);
    b=x3(p,3);
end
if bdp==0.25
    c=x2(p,2);
    b=x2(p,3);
end
