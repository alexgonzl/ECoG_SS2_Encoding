function out = getBinaryPerf(y,yhat,posClass,negClass)

out.TP      = find(y==posClass & yhat==posClass);
out.nTP     = numel(out.TP);

out.TN      = find(y==negClass & yhat==negClass);
out.nTN     = numel(out.TN);

out.FP      = find(y==negClass & yhat==posClass);
out.nFP     = numel(out.FP);

out.FN      = find(y==posClass & yhat==negClass);
out.nFN     = numel(out.FN);

out.AC     = (out.nTP+out.nTN)/(out.nTP+out.nTN+out.nFP+out.nFN);

out.SE      = out.nTP/(out.nTP+out.nFN);
out.SP      = out.nTN/(out.nTN+out.nFP);