function  plottingWrapper(figNum)
% plots for ECoG_SS2enconding. 

addpath ~/Documents/ECOG_Mem_Encoding/lib/
addpath ~/Documents/ECOG_Mem_Encoding/Analysis/

switch figNum
	case 1
	% I. Behavioral
	behavPlots();
	case 2
	% II. Activity Plots
	% ROI Encoding Activity 
	SS2e_MBAnalysis('preStim2',1)

	otherwise
		error('figure has not been implemented.')
end
