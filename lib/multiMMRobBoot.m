
function bootres = multiMMRobBoot(Y,Z,Beta_0,Sigma_0,MMBeta, MMSigma, R, bdp)

% robust bootstrap for multivariate MM regression or location/shape estimation 
% INPUT:
%   Y : n x q response matrix
%   Z : n x p covariates matrix (ones(n,1) for location/shape estimation)
%   MMBeta : p x q MM-coefficients
%   MMSigma : q x q MM-covariance
%   Beta_0 : p x q S-coefficients
%   Sigma_0 : q x q S-covariance
%   R : number of bootstrap samples
%   bdp : breakdown point of S-estimate (probably =0.5)
%
% OUTPUT: 
%   bootres.centered: (2*dimens x R) centered recomputations of MM and S-estimates:
%                             - first p*q rows: MM location (or regression)      
%                             - next q*q rows : MM shape matrix 
%                             - next q*q rows : S covariance matrix
%                             - final p*q rows: S location (or regression)
%                   (all in vec-form, columns stacked on top of each other) 
%   
%   bootres.eigGcentered: (q x R) centered recomputations of EIGENVALUES of
%                                   MM shape estimate
%   bootres.eigveccentered: (q*q x R) centered recomputations of 
%                           EIGENVECTORS of MM shape estimate (in vec-form)
%
%   bootres.pvarcentered: (q-1 x R) centered recomputations of PERCENTAGES
%                                         of variance for MM shape estimate
%
%   bootres.angles: (q x R) angles that bootstrap eigenvectors have with
%                                   original MM-eigenvectors (in radians)
%
%   bootres.MMest: (2*dimens x 1) original MM (and S) estimates in vec-form

[n q]=size(Y);
[n p]=size(Z);
dimens=p*q+q*q;

Iq=diag(ones(q,1));
Ip=diag(ones(p,1));

c0 = Tbsc(bdp,q);
b = (c0/6) * Tbsb(c0,q);
c1 = Tbsc1(.95,q);

S_0inv=inv(Sigma_0);  % nodig?
auxscalesq=det(Sigma_0)^(1/q);
auxscale=sqrt(auxscalesq);
MMGamma=auxscalesq^(-1)*MMSigma;
MMGinv=inv(MMGamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        calculate jacobian                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%           p*q          q*q          q*q            p*q  
%       -----------------------------------------------------
%       |           |            |            |             |
%  p*q  | g1_MMBeta | g1_MMGamma | g1_Sigma_0 |      0      |
%       |           |            |            |             |
%       -----------------------------------------------------
%       |           |            |            |             |
%  q*q  | g2_MMBeta | g2_MMGamma | g2_Sigma_0 |      0      |
%       |           |            |            |             |  
%       -----------------------------------------------------
%       |           |            |            |             |
%  q*q  |     0     |      0     | g3_Sigma_0 |  g3_Beta_0  |
%       |           |            |            |             |
%       -----------------------------------------------------
%       |           |            |            |             |
%  p*q  |     0     |      0     | g4_Sigma_0 |  g4_Beta_0  |
%       |           |            |            |             |
%       -----------------------------------------------------
%             1            2            3            4 

% first fill up lower right part: the S-part : g3, g4

restildematrix=Y-Z*Beta_0;
ditildevec=sqrt(mahsq(restildematrix, zeros(1,q), Sigma_0, n));
ditildevec(ditildevec(:)<1e-5)=1e-5;
uditildevec=rhobiweightder1(ditildevec,c0)./ditildevec;
wditildevec=(rhobiweightder2(ditildevec,c0).*ditildevec-rhobiweightder1(ditildevec,c0))./ditildevec.^3;
zditildevec=rhobiweightder2(ditildevec,c0);
wwditildevec=rhobiweightder1(ditildevec,c0).*ditildevec-rhobiweight(ditildevec,c0);

uditildemat=diag(uditildevec);
Atilde=Z'*uditildemat*Z;
Btilde=Z'*uditildemat*Y;
% Vtilde=q*restildematrix'*uditildemat*restildematrix;

term1a=zeros(p*p,p*q);
term1b=zeros(p*q,p*q);
term2a=zeros(q*q,p*q);
term2b=zeros(q*q,p*q);
term2c=zeros(1,p*q);
term3a=zeros(p*p,q*q);
term3b=zeros(p*q,q*q);
term4a=zeros(q*q,q*q);
term4b=zeros(1,q*q);

for i=1:n
    Zi=Z(i,:)';
    Yi=Y(i,:)';
    resi=restildematrix(i,:)';
    vecZiZi=vecop(Zi*Zi');
    vecZiYi=vecop(Zi*Yi');
    vecresiresi=vecop(resi*resi');
    wdi=wditildevec(i);
    zdi=zditildevec(i);
    udi=uditildevec(i);
    tveci=vecop(Zi*resi'*S_0inv)';
    tvecSi=vecop(S_0inv*resi*resi'*S_0inv)';

    term1a=term1a+wdi*(vecZiZi*tveci);
    term1b=term1b+wdi*(vecZiYi*tveci);

    term2a=term2a+wdi*(vecresiresi*tveci);
    term2b=term2b+udi*((kron(Iq,resi)+kron(resi,Iq))*(kron(Zi',Iq)*commut(p,q))); 
    term2c=term2c+zdi*tveci;
    
    term3a=term3a+wdi*(vecZiZi*tvecSi);
    term3b=term3b+wdi*(vecZiYi*tvecSi);

    term4a=term4a+wdi*(vecresiresi*tvecSi);
    term4b=term4b+zdi*tvecSi;
end

Atildeinv=inv(Atilde);

partder1=(kron(Btilde,Ip)'*kron(Atildeinv',Atildeinv)*term1a)-(kron(Iq,Atildeinv)*term1b);
partder2=-q/(b*n)*(term2a+term2b)+1/(b*n)*vecop(Sigma_0)*term2c;
partder3=1/2*(kron(Btilde,Ip)'*kron(Atildeinv',Atildeinv)*term3a)-1/2*(kron(Iq,Atildeinv)*term3b);
partder4=-q/(2*b*n)*term4a+1/(2*b*n)*vecop(Sigma_0)*term4b-1/(b*n)*sum(wwditildevec)*diag(ones(q*q,1));

Part3_3=partder4;
Part3_4=partder2;
Part4_3=partder3;
Part4_4=partder1;

% end S-part

% now g1, g2 

resmatrix=Y-Z*MMBeta;
divec=sqrt(mahsq(resmatrix, zeros(1,q), MMGamma, n));
divec(divec(:)<1e-5)=1e-5;
udivec=rhobiweightder1(divec/auxscale,c1)./divec;
wdivec=(rhobiweightder2(divec/auxscale,c1).*divec/auxscale-rhobiweightder1(divec/auxscale,c1))./divec.^3;
vdivec=rhobiweightder2(divec/auxscale,c1);

udimat=diag(udivec);
A=Z'*udimat*Z;
B=Z'*udimat*Y;
V=resmatrix'*udimat*resmatrix;

termg1_MMBeta_a=zeros(p*p,p*q);
termg1_MMBeta_b=zeros(p*q,p*q);
termg1_MMGamma_a=zeros(p*p,q*q);
termg1_MMGamma_b=zeros(p*q,q*q);
termg2_MMBeta_a=zeros(q*q,p*q);
termg2_MMBeta_b=zeros(q*q,p*q);
termg2_MMGamma=zeros(q*q,q*q);
termg1_MMSigma_0a=zeros(p*p,1);
termg1_MMSigma_0b=zeros(p*q,1);
termg2_MMSigma_0=zeros(q*q,1);

for i=1:n
    Zi=Z(i,:)';
    Yi=Y(i,:)';
    resi=resmatrix(i,:)';
    vecZiZi=vecop(Zi*Zi');
    vecZiYi=vecop(Zi*Yi');
    vecresiresi=vecop(resi*resi');
    wdi=wdivec(i);
    udi=udivec(i);
    vdi=vdivec(i);
    tveci=vecop(Zi*resi'*MMGinv)';
    tvecSi=vecop(MMGinv*resi*resi'*MMGinv)';

    termg1_MMBeta_a=termg1_MMBeta_a+wdi*(vecZiZi*tveci);
    termg1_MMBeta_b=termg1_MMBeta_b+wdi*(vecZiYi*tveci);

    termg1_MMGamma_a=termg1_MMGamma_a+wdi*(vecZiZi*tvecSi);
    termg1_MMGamma_b=termg1_MMGamma_b+wdi*(vecZiYi*tvecSi);

    termg2_MMBeta_a=termg2_MMBeta_a+wdi*(vecresiresi*tveci);
    termg2_MMBeta_b=termg2_MMBeta_b+udi*((kron(Iq,resi)+kron(resi,Iq))*(kron(Zi',Iq)*commut(p,q))); 

    termg2_MMGamma=termg2_MMGamma+wdi*(vecresiresi*tvecSi);
    
    termg1_MMSigma_0a=termg1_MMSigma_0a+vdi*vecZiZi;
    termg1_MMSigma_0b=termg1_MMSigma_0b+vdi*vecZiYi;

    termg2_MMSigma_0=termg2_MMSigma_0+vdi*vecresiresi;
    
end

Ainv=inv(A);

Part1_1=(kron(B,Ip)'*kron(Ainv',Ainv)*termg1_MMBeta_a)-(kron(Iq,Ainv)*termg1_MMBeta_b);
Part1_2=1/2*(kron(B,Ip)'*kron(Ainv',Ainv)*termg1_MMGamma_a)-1/2*(kron(Iq,Ainv)*termg1_MMGamma_b);
Part2_1=-det(V)^(-1/q)*(diag(ones(q*q,1))-1/q*vecop(V)*vecop(inv(V)')')*(termg2_MMBeta_a+termg2_MMBeta_b);
Part2_2=-1/2*det(V)^(-1/q)*(diag(ones(q*q,1))-1/q*vecop(V)*vecop(inv(V)')')*termg2_MMGamma;
Part1_3=-1/2/q/auxscale*(kron(B,Ip)'*kron(Ainv',Ainv)*termg1_MMSigma_0a-kron(Iq,Ainv)*termg1_MMSigma_0b)*vecop(S_0inv')';
Part2_3=-1/2/q/auxscale*det(V)^(-1/q)*(diag(ones(q*q,1))-1/q*vecop(V)*vecop(inv(V)')')*termg2_MMSigma_0*vecop(S_0inv')';

Part1_4=zeros(p*q,p*q);
Part2_4=zeros(q*q,p*q);
Part3_1=zeros(q*q,p*q);
Part3_2=zeros(q*q,q*q);
Part4_1=zeros(p*q,p*q);
Part4_2=zeros(p*q,q*q);

jacobian=[Part1_1 Part1_2 Part1_3 Part1_4;
            Part2_1 Part2_2 Part2_3 Part2_4; 
            Part3_1 Part3_2 Part3_3 Part3_4;
            Part4_1 Part4_2 Part4_3 Part4_4];

Idim=diag(ones(dimens*2,1));
lincorrectmat=inv(Idim-jacobian);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% put all estimates (coefs and covariances) in one column 
vecestim=zeros(dimens*2,1);
vecestim(1:(p*q))=vecop(MMBeta);
vecestim((p*q+1):dimens)=vecop(MMGamma);
vecestim((dimens+1):(dimens+(q*q)))=vecop(Sigma_0);
vecestim((dimens+q*q+1):(dimens*2))=vecop(Beta_0);

% compute corresponding original eigenvalue/eigenvector estimates
[VGamma DGamma]=eig(MMGamma);
[VSigma DSigma]=eig(MMSigma);
[seigMMGamma IXGamma]=sort(diag(DGamma));
[seigMMSigma IXSigma]=sort(diag(DSigma));
eigenvectorsMM=VSigma(:,IXSigma);
majorloading=zeros(1,q);
for k=1:q
    [dummy IX]=sort(abs(eigenvectorsMM(:,k)));
    majorloading(k)=IX(q);
    eigenvectorsMM(:,k)=sign(eigenvectorsMM(IX(q),k))*eigenvectorsMM(:,k);
end    
varperc=zeros(q-1,1);
for k=1:(q-1)
    varperc(k)=sum(seigMMGamma((q-k+1):q))/sum(seigMMGamma);
end

% to draw bootstrap samples
[dummy bootmatrix]=bootstrp(R,'mean',zeros(n,1));

bootbiasmat = zeros(dimens*2,R);  
booteigvalsGamma = zeros(q,R);
booteigvecs = zeros(q*q,R);
bootpercvars = zeros(q-1,R);
bootangles = zeros(q,R);
bootsampleOK = ones(1,R);

for r=1:R
    Yst=Y(bootmatrix(:,r),:);
    Zst=Z(bootmatrix(:,r),:);
    resmatrixst=Yst-Zst*MMBeta;
    divecst=sqrt(mahsq(resmatrixst,zeros(1,q),MMSigma,n));
    divecst(divecst(:)<1e-5)=1e-5;
    udivecst=rhobiweightder1(divecst,c1)./divecst;
    
    restildematrixst=Yst-Zst*Beta_0;
    ditildevecst=sqrt(mahsq(restildematrixst,zeros(1,q),Sigma_0,n));
    ditildevecst(ditildevecst(:)<1e-5)=1e-5;
    uditildevecst=rhobiweightder1(ditildevecst,c0)./ditildevecst;
    wwditildevecst=rhobiweightder1(ditildevecst,c0).*ditildevecst-rhobiweight(ditildevecst,c0);   
    
    udimatst=diag(udivecst);
    Bst=inv(Zst'*udimatst*Zst)*(Zst'*udimatst*Yst);
    Gst=resmatrixst'*udimatst*resmatrixst;
    Gst=det(Gst)^(-1/q)*Gst;    
%     Vst=auxscalesq*Gst;
    
    uditildematst=diag(uditildevecst);
    V0st_term1=1/(b*n)*q*restildematrixst'*uditildematst*restildematrixst;
    V0st_term2=1/(b*n)*sum(wwditildevecst)*Sigma_0;
    V0st=V0st_term1-V0st_term2;
    
    B0st=inv(Zst'*uditildematst*Zst)*(Zst'*uditildematst*Yst);
    
    % list uncorrected bootstrap recomputations
    vecfst=zeros(dimens*2,1);
    vecfst(1:(p*q))=vecop(Bst);
    vecfst((p*q+1):dimens)=vecop(Gst);
    vecfst((dimens+1):(dimens+q*q))=vecop(V0st);
    vecfst((dimens+q*q+1):(dimens*2))=vecop(B0st);
    
    % compute centered, corrected fast bootstrap estimates
    fstbias=vecfst-vecestim;
    bootbiasmat(:,r)=lincorrectmat*fstbias;  
 
    % compute bootstrap eigenvalues/vectors, percentages and angles
    correctedMMGammast=reconvec(bootbiasmat((p*q+1):dimens,r)+vecestim((p*q+1):dimens),q);
    [VGammast DGammast]=eig(correctedMMGammast);
    [eigenvaluesGammast IXGammast]=sort(diag(DGammast));
    eigenvectorsGammast=VGammast(:,IXGammast);
    for k=1:q
        eigenvectorsGammast(:,k)=sign(eigenvectorsGammast(majorloading(k),k))*eigenvectorsGammast(:,k);
    end    
    if any(eigenvaluesGammast<0)
        bootsampleOK(r)=0;
    end
    
    percvarGammast=zeros(q-1,1);
    for k=1:(q-1)
        percvarGammast(k)=sum(eigenvaluesGammast((q-k+1):q))/sum(eigenvaluesGammast);
    end

    Svecs=eigenvectorsGammast;
    for k=1:q
        bootangles(k,r)=acos(abs(Svecs(:,k)'*eigenvectorsMM(:,k)/sqrt(Svecs(:,k)'*Svecs(:,k))/sqrt(eigenvectorsMM(:,k)'*eigenvectorsMM(:,k))));
    end
    booteigvalsGamma(:,r)=eigenvaluesGammast-seigMMGamma;
    booteigvecs(:,r)=vecop(eigenvectorsGammast)-vecop(eigenvectorsMM);
    bootpercvars(:,r)=percvarGammast-varperc;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


bootindices = 1:R;
bootindicesOK = bootindices(bootsampleOK==1);

bootres.centered=bootbiasmat(:,bootindicesOK);
bootres.eigGcentered=booteigvalsGamma(:,bootindicesOK);
bootres.eigveccentered=booteigvecs(:,bootindicesOK);
bootres.pvarcentered=bootpercvars(:,bootindicesOK);
bootres.angles=bootangles(:,bootindicesOK);
bootres.MMest=vecestim;

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

function recon=reconvec(vector,ncol)
% reconstructs vecop'd matrix
lcol=length(vector)/ncol;
rec=zeros(lcol,ncol);
for i=1:ncol
    rec(:,i)=vector(((i-1)*lcol+1):(i*lcol));
end
recon=rec;

%-------------------------------------------------------------------------- 

function res=gint(k,c,p)

% this procedures computes the integral from zero to c of
% the function r^k g(r^2), where g(||x||^2) is the density function
% of a p-dimensional standardnormal distribution 
% 
% Direct calls: none
e=(k-p-1)/2;
numerator=(2^e)*gamcdf((c^2)/2,(k+1)/2)*gamma((k+1)/2);
res=(numerator/(pi^(p/2)));

%-------------------------------------------------------------------------- 

function res=Tbsb(c,p)
% Direct calls: gint
y1=gint(p+1,c,p)/2-gint(p+3,c,p)/(2*c^2)+gint(p+5,c,p)/(6*c^4);
y2=(6/c)*2*(pi^(p/2))/gamma(p/2);
y3=c*(1-chi2cdf(c^2,p));
res=(y1*y2+y3);

%-------------------------------------------------------------------------- 

function res = Tbsc(alpha,p)
% constant for Tukey Biweight S 
%
% Direct cals: Tbsb
talpha = sqrt(chi2inv(1-alpha,p));
maxit = 1000; 
eps = 10^(-8);
diff = 10^6;
ctest = talpha;
iter = 1;
while ((diff>eps) && iter<maxit)
    cold = ctest;
    ctest = Tbsb(cold,p)/alpha;
    diff = abs(cold-ctest);
    iter = iter+1;
end
res = (ctest);

%-------------------------------------------------------------------------- 

function res=sigma1(c,p)
% Direct calls: gint

Cp = 2*(pi^(p/2))/gamma(p/2);
gamma1_1 = gint(p+1,c,p) - 6*gint(p+3,c,p)/(c^2) + 5*gint(p+5,c,p)/(c^4);
gamma1_2 = gint(p+1,c,p) - 2*gint(p+3,c,p)/(c^2) + gint(p+5,c,p)/(c^4);
gamma1 = Cp * ( gamma1_1 + (p+1)*gamma1_2 ) / (p+2);

s1 = Cp * ( gint(p+3,c,p) - 4*gint(p+5,c,p)/(c^2) + 6*gint(p+7,c,p)/(c^4) - 4*gint(p+9,c,p)/(c^6) + gint(p+11,c,p)/(c^8) ) / gamma1^2 * p/(p+2);

res=s1;

%-------------------------------------------------------------------------- 

function res = Tbsc1(shape_eff,p)
% constant for second Tukey Biweight rho-function for MM, for fixed shape-efficiency 
%
% Direct cals: sigma1, Tbsc

maxit = 1000; 
eps = 10^(-8);
diff = 10^6;
ctest = Tbsc(.5,p);
iter = 1;
while ((diff>eps) && iter<maxit)
    cold = ctest;
    ctest = cold*shape_eff*sigma1(cold,p);
    diff = abs(cold-ctest);
    iter = iter+1;
end
res = (ctest);

%-------------------------------------------------------------------------- 

function Kpm=commut(p,m)

%computes commutation matrix
%p = no of rows
%m = no of columns (of matrix which follows Kpm)

kompm = zeros(p*m,p*m);
for k=1:(p*m)
    l=(k-1-(ceil(k/m)-1)*m)*p+ceil(k/m);
    kompm(k,l)=1;
end

Kpm = kompm;

    
    





