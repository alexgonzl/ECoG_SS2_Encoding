
function res=multiMM(Y)

% computes multivariate location and shape (M)M-estimates with auxiliary S-scale 
% input:
%   Y = data
% output:
%   res.loc : MM-location estimate
%   res.shape : MM-shape matrix
%   res.covariance : MM-covariance
%   res.Sloc : S-location estimate
%   res.Sshape : S-shape matrix
%   res.auxscale : S-scale
% breakdown point (bdp) set at .5
% efficiciency : 95% shape-efficiency
%
% (this file contains its own version of the fastS-code: 'multiS')

[n q] = size(Y);
c = Tbsc1(.95,q);

Snsamp = 100; % number of elemental starts in S-estimation

% first compute 50% breakdown S-estimator
Sresult = multiS(ones(n,1), Y, Snsamp, 0.5);

auxscale = Sresult.scale;
newG = Sresult.Gamma;
newB = Sresult.Beta;
newR = Y-ones(n,1)*newB;
psres = sqrt(mahsq(newR, zeros(1,q), newG, n));
newobj = mean(rhobiweight(psres/auxscale,c));
origobj = newobj;

maxiter = 20; % maximum number of IRLS steps
mtol = 1e-7;

% compute M-estimate with auxilliary scale through IRLS steps, starting
% from S-estimate
iteration=1;
oldobj = newobj + 1;
while ((oldobj - newobj) > mtol) && (iteration < maxiter)
    oldobj = newobj;
    w = psibiweight(psres/auxscale,c)./psres;
    wbig = diag(w);
    newB = sum(wbig*Y)/sum(w);
    newG = q*(newR'*wbig*newR); 
    newG = det(newG)^(-1/q)*newG;
    newR = Y-ones(n,1)*newB;
    psres = sqrt(mahsq(newR, zeros(1,q), newG, n));
    newobj = mean(rhobiweight(psres/auxscale,c));
    iteration = iteration+1;
end

if newobj <= origobj
    res.loc = newB;
    res.shape = newG;
    res.covariance = newG*auxscale^2;
else % isn't supposed to happen
   res.loc = Sresult.Beta;
   res.shape = Sresult.Gamma;
   res.covariance = Sresult.Gamma*auxscale^2;
end

res.Sloc = Sresult.Beta;   
res.Sshape = Sresult.Gamma;
res.auxscale = auxscale;

%--------------------------------------------------------------------------

function result = multiS(X,Y,nsamp,bdp)

% fast S algorithm for multivariate regression or location estimation
% Y = response matrix (n x m)
% X = covariates matrix (n x p), possibly including intercept column
% (X = ones(n,1) in case of location estimation)
% nsamp = number of elemental starts, e.g. 20 
% bdp = breakdown point (0.5 or 0.25 or 0.15)

seed=0;
[n p]=size(X);
[n m]=size(Y);

loop=1;
c = Tbsc(bdp,m);
b = (c/6) * Tbsb(c,m);
bestr= 5; % number of best solutions to keep for further C-steps
k = 3; % number of C-steps on elemental starts
bestbetas = zeros(p*m, bestr);
bestgammas = zeros(m*m, bestr);
bestscales = 1e20 * ones(bestr,1);
sworst = 1e20;

while loop<=nsamp
    % find a (p+m)-subset in general position.
    rankxy=0;
    while rankxy<(p+m)
        [ranset,seed]=randomset(n,p+m,seed);
        Xj=X(ranset,:);
        Yj=Y(ranset,:);
        rankxy=rank([Xj Yj]);
    end
    Bj=Xj\Yj;
    Rj=Yj-Xj*Bj;
    Sj=Rj'*Rj/(p+m-1);
    Gj=det(Sj)^(-1/m)*Sj;
    % perform k steps of IRLS on elemental start
    res=IRLSstep(X, Y, Bj, Gj, 0, k, c, b);
    Betarw = res.Beta;
    Gammarw = res.Gamma;
    scalerw = res.scale;
    psresrw = sqrt(mahsq(Y-X*Betarw,zeros(1,m),Gammarw,n));
    if loop > 1
        % check whether new Beta and new Gamma belong to the top best Betas; if so keep
        % Beta and Gamma with corresponding scale.
        if mean(rhobiweight(psresrw/sworst,c)) < b
            [yss,yi]=sort(bestscales);
            ind=yi(bestr);
            bestscales(ind) = scale1(psresrw,b,c,scalerw);
            bestbetas(:,ind) = vecop(Betarw);
            bestgammas(:,ind) = vecop(Gammarw);
            sworst = max(bestscales);
        end
    else
        bestscales(bestr) = scale1(psresrw,b,c,scalerw);
        bestbetas(:,bestr) = vecop(Betarw);
        bestgammas(:,bestr) = vecop(Gammarw);
    end
    loop=loop+1;
end

[superbestscale ibest] = max(bestscales);
superbestbeta = reconvec(bestbetas(:,ibest),m);
superbestgamma = reconvec(bestgammas(:,ibest),m);

% perform C-steps on best 'bestr' solutions, until convergence (or maximum 50 steps) 
for i=bestr:-1:1
    tmp = IRLSstep(X,Y,reconvec(bestbetas(:,i),m), reconvec(bestgammas(:,i),m), bestscales(i), 50, c, b);
    if tmp.scale < superbestscale
        superbestscale = tmp.scale;
        superbestbeta = tmp.Beta;
        superbestgamma = tmp.Gamma;
    end
end

result.Beta=superbestbeta;
result.Gamma=superbestgamma;
result.scale=superbestscale;


% -------------------------------------------------------------------

function res=IRLSstep(X,Y,initialBeta, initialGamma, initialscale, k, c, b)
  
convTol = 1e-10;

[n,m]=size(Y); 

Beta=initialBeta;
Res=Y-X*Beta;
psres = sqrt(mahsq(Res,zeros(1,m),initialGamma,n));
if initialscale > 0
    scale = initialscale;
else
    scale = median(psres)/.6745;
end

iter = 0;
betadiff = 1;

while ( (betadiff > convTol) && (iter < k))
    iter = iter + 1;
    scale = sqrt(scale^2 * mean(rhobiweight(psres/scale,c))/b);
    w = scaledpsibiweight(psres/scale,c);
    wbig = diag(w);
    newBeta = (X'*wbig*X)\(X'*wbig*Y);
    newGamma = m*(Res'*wbig*Res);  % computing Gamma with newBeta already? (instead of with Beta)
    newGamma = det(newGamma)^(-1/m)*newGamma;
    Res = Y-X*newBeta;
    betadiff = norm(newBeta-Beta,1)/norm(Beta,1);
    Beta = newBeta;
    psres = sqrt(mahsq(Res,zeros(1,m),newGamma,n));
end

res.Beta = newBeta;
res.Gamma = newGamma;  
res.scale = scale;

%--------------------------------------------------------------------------  

function sc = scale1(u, b, c, initialsc) 
% from Kristel's fastSreg
if nargin<3
    initialsc = median(abs(u))/.6745;
end
max.it = 100;
sc = initialsc;
i = 0; 
eps = 1e-10;
err = 1;
while  (( i < max.it ) && (err > eps))
    sc2 = sqrt( sc^2 * mean(rhobiweight(u/sc,c)) / b);
    err =abs(sc2/sc - 1);
    sc = sc2;
    i=i+1;
end

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

% --------------------------------------------------------------------

function rho=rhobiweight(x,c)

% Computes Tukey's biweight rho function met constante c voor alle waarden
% in de vector x.

hulp=x.^2/2-x.^4/(2*c^2)+x.^6/(6*c^4);
rho=hulp.*(abs(x)<c)+c^2/6.*(abs(x)>=c);

% --------------------------------------------------------------------

function psi=psibiweight(x,c)

% Computes Tukey's biweight psi function met constante c voor alle waarden
% in de vector x.

hulp=x-2.*x.^3/(c^2)+x.^5/(c^4);
psi=hulp.*(abs(x)<c);

% --------------------------------------------------------------------

function psi=scaledpsibiweight(x,c)

% Computes Tukey's biweight psi function divided elementwise with x,
% with constant c for all values in the vector x.
hulp=1-2.*x.^2/(c^2)+x.^4/(c^4);
psi=hulp.*(abs(x)<c);

% --------------------------------------------------------------------

function [random,seed]=uniran(seed)

% The random generator.

seed=floor(seed*5761)+999;
quot=floor(seed/65536);
seed=floor(seed)-floor(quot*65536);
random=seed/65536.D0;

% --------------------------------------------------------------------


function [ranset,seed] = randomset(tot,nel,seed)

for j = 1:nel
   [random,seed]=uniran(seed);
   num=floor(random*tot)+1;
   if j > 1
      while any(ranset==num)
         [random,seed]=uniran(seed);
         num=floor(random*tot)+1;
      end
   end
   ranset(j)=num;
end

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


