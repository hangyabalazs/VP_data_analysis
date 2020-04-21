function acg_sorting(cellids, sourcedir, resdir_burst, resdir_nonburst)
%ACG_SORTING(CELLIDS, SOURCEDIR, RESDIR_BURST, RESDIR_NONSCBURST) finds bursting neurons based  on their autocorrelograms.
%
%ACG_SORTING(CELLIDS, SOURCEDIR, RESDIR_BURST, RESDIR_NONSCBURST) finds
%bursting neurons based on their autocorrelograms. Burst index (BI) was 
%calculated by the normalized difference between maximum ACG for lags  
%0-10 ms and mean ACG for lags 180-200 ms, where the normalizing factor was
%the greater of the two numbers, yielding an index between -1 and 1 (Royer et al., 2012)
%Neurons with Burst Index >0.2 are considered as bursting neurons and
%their autocorrelogram is saved to RESDIR_BURST. Acg of non-bursting
%neurons are saved to RESDIR_NONBURST.

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   14-Apr-2020

if ~isfolder(resdir_burst)
    mkdir(resdir_burst)
end

if ~isfolder(resdir_nonburst)
    mkdir(resdir_nonburst)
end

numPC = length(cellids);
load([sourcedir '\ACG_matrices_.mat']); %load acg matrix 

for k = 1:numPC
    cellid = char(cellids{k});
    sttc = regexprep(cellid,'\.','_');  % string to compare
    
    % Move acg to folder
    cd(sourcedir);
    if BurstIndex(k) > 0.2  % copy bursting cell acgs to folder
        copyfile(['*' sttc '*'], resdir_burst);
    elseif ~isnan(BurstIndex(k)) && BurstIndex(k)<=0.2
        copyfile(['*' sttc '*'], resdir_nonburst);
    end
end