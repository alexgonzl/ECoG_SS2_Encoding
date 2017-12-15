function out = cellTStatWrapper(X)

N = numel(X);
out.uniVariateT = zeros(N,1);
out.uniVariateP = zeros(N,1);
for ii =1:N
    [~,p,~,t] = ttest(X{ii});
    out.uniVariateT(ii) = t.tstat;
    out.uniVariateP(ii) = p;
end

out.nCombs      = nchoosek(N,2);
out.Combs       = zeros(out.nCombs,2);
out.biVariateT  = zeros(out.nCombs,1);
out.biVariateP  = zeros(out.nCombs,1);
cnt = 1;
for ii = 1:N
    for jj = (ii+1):N
        out.Combs(cnt,:)       = [ii jj];
        [~,p,~,t] = ttest2(X{ii},X{jj});
        out.biVariateT(cnt)  = t.tstat;
        out.biVariateP(cnt)  = p;
        cnt = cnt +1;
    end
end


p = [out.uniVariateP; out.biVariateP];
[~,a] = fdr_bh(p);

out.FDR_Thr = a;

