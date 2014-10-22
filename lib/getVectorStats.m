function out = getVectorStats(data)

data(isnan(data))=[];

measures = {'numel','mean','median','se','std','kurtosis'};

for mm = 1:numel(measures)
    out.(measures{mm}) = eval([measures{mm} '(data)']);
end

out.quant25     = quantile(data,0.25);
out.quant75     = quantile(data,0.75);



