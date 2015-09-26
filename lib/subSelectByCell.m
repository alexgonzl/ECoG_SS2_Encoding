
function out = subSelectByCell(cellMat,cellSelection)
% utility function for selecting trials in a cell structure 
% that holds ecog data for multiple subjects, channels, trials, samples

nCells      = numel(cellMat);
dims        = max(cellfun('ndims',cellMat));
cellSizes   = zeros(dims,nCells);

for dd = 1:dims
    cellSizes(dd,:) = cellfun('size',cellMat,dd);
end

out = cell(nCells,1);
for ii = 1:nCells
	for jj = 1:cellSizes(1,ii)
        if sum(cellSelection{ii})>1  
            out{ii}(jj,:,:) = cellMat{ii}(jj,cellSelection{ii},:);
        elseif sum(cellSelection{ii})==1
            out{ii}(jj,1,:) = cellMat{ii}(jj,cellSelection{ii},:);
        else
            out{ii}(jj,1,:) = nan;
        end
	end
end