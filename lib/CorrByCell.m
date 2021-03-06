

function out = CorrByCell(cellMat,cellRTs)
% utility function

nCells      = numel(cellMat);
dims        = max(cellfun('ndims',cellMat));
cellSizes   = zeros(dims,nCells);

for dd = 1:dims
    cellSizes(dd,:) = cellfun('size',cellMat,dd);
end

N = sum(cellSizes(1,:));

out = cell(N,1);
cnt = 1;
for ii = 1:nCells    
   for jj = 1:cellSizes(1,ii)
       X= squeeze(cellMat{ii}(jj,:,:));
       out{cnt} = corr(X,cellRTs{ii},'rows','complete','type','spearman');
       cnt = cnt +1;
   end
end