function [phasic_cells, BetaIndex, GammaIndex] = pr_freq_detection_VP(vpcells, resdir, resdir2, resdir3, acgsource, psthsource)
%PR_FREQ_DETECTION_VP finds rhythmic neurons based  on their autocorrelograms.
%
%   [PHASIC_CELLS, BETAINDEX, GAMMAINDEX] = PR_FREQ_DETECTION_VP(RESDIR,
%   ACGSOURCE, PSTHSOURCE) finds rhythmic neurons within the beta (16-30
%   Hz) and gamma (30-100Hz) frequency rage based on their autocorrelogram
%   peaks. Bursting neurons (acg peak with <10 ms lag) and slow rhytmic
%   cells (>60 ms lag) are excluded from the analysis. To determine
%   rhythmicity, PR_FREQ_DETECTION uses the autocorrelation matrix made by
%   VPACG (ACG_matrices.mat). PR_FREQ_DETECTION detects peaks on the acg in
%   the defined (beta and gamma) range and calculates a ratio between the
%   peak vs. baseline (BetaIndex, GammaIndex).
%
%   Rhythmic cell autocorrelograms and response PSTHs are selected from the
%   PSTHs (PSTHSOURCE) and autocorrelograms (ACGSOURCE) made by VPACG and
%   VPRESPONSESORTER respectively.
%
%   See also VPACG and VPRESPONSESORTER.

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   05-Feb-2020 

%   Code review: BH 4/8/20

% Directories
if ~isfolder(resdir)
    mkdir(resdir);
end
if ~isfolder(resdir2)
    mkdir(resdir2);
end
if ~isfolder(resdir3)
    mkdir(resdir3);
end

% Load ACG matrices - necessary for rhythmicity detection
load(fullfile(acgsource,'ACG_matrices.mat')); %#ok<*LOAD>

% Select cells
if isempty(vpcells)
    vpcells = cellids; %if no cellid inputs - use the ones in the matrix
end

% Preallocate
numCells = length(vpcells);
[peaks, peak_inx, lagms] = deal(nan(1,numCells));

% Exclude bursting neurons
res = 0.5;  % ACG resolution
phasewin = 10;
ACGh = CCR(:,1001:2000);  % half acg
lagsh = lags(1001:2000); % corresponding half of time vector in ms
for i = 1:length(cellids)  % loop through cells
    cACGh = ACGh(i,:);
    peaks(i) = max(cACGh);
    peak_inxx = find(cACGh == peaks(i));
    peak_inx(i)= peak_inxx(end);
    lagms(i) = lagsh(peak_inx(i));
end
phasic_cells = cellids(10<lagms & lagms<=60); % only choosing beta and gamma cells and eliminate bursting ones (lag>10ms)

% Beta and gamma index calculation
numPC = length(phasic_cells); % number of phasic cells
[betainx, gammainx, cellinx] = deal(nan(1,numPC));
for j = 1:numPC  % loop through non-bursting cells with <60 ms ACG peaks
    cellid = phasic_cells{j};
    cellinx(j) = find(strcmp(cellids,cellid)); % getting index from original cellidlist
    pACGh = ACGh(cellinx(j),:);  % acg of the current phasic cell
    
    % Find beta peak
    b_peak = pACGh(lagsh<=60 & lagsh>=30);
    [b_peak_val, binx] = max(b_peak);
    flank = phasewin / res;
    mbp = mean(pACGh(59+binx-flank:59+binx+flank)); % mean beta peak
    
    % Calculate baseline and Beta Index
    peak_loc = 59+binx;
    bsl = (pACGh(round(peak_loc/2))+pACGh(peak_loc*2)) / 2;  %baseline
    betainx(j) = (mbp-bsl)/max(mbp,bsl); % beta index
    
    % Find gamma peak
    g_peak = pACGh(lagsh<=31 & lagsh>=10);
    g_peak_inx = find(lagsh<=31 & lagsh>=10);
    [g_peak_val, ginx] = max(g_peak);
    flank = phasewin / res;
    mgp = mean(pACGh(20+ginx-flank:20+ginx+flank)); % mean gamma peak
    
    % Calculate baseline and Gamma Index
    peak_loc2 = 20 + ginx;
    bsl2 = (pACGh(round(peak_loc2/2))+pACGh(peak_loc2*2)) / 2;
    gammainx(j) = (mgp-bsl2)/max(mgp,bsl2); % gamma index
end
BetaIndex = betainx;
GammaIndex = gammainx;
ph_Refractory = Refractory(cellinx);   % refractory period

% Copy psth and acg plots of rhythmic cells
filenames1 = sourcefiles(acgsource); % acg filenames
filenames2 = sourcefiles(psthsource); % psth filenames
for k = 1:numPC
    cellid = char(phasic_cells{k});
    sttc = regexprep(cellid,'\.','_');  % string to compare
    
    % Move acg to folder
    cd(acgsource);
    if BetaIndex(k) > 0.4 || GammaIndex(k) > 0.25  % copy rhythmic cell acgs to folder
        copyfile(['*' sttc '*'], resdir);
    else
        copyfile(['*' sttc '*'], resdir2);  % copy non rhythmic cell acgs to another folder
    end

    % Move psth to folder
    cd(psthsource);
    if BetaIndex(k) > 0.4 || GammaIndex(k) > 0.25  % copy rhythmic cell acgs to folder
        copyfile(['*' sttc '*'], resdir3);
    end
end

% Save
fnmm = 'phasic_response.mat'; % save phasic cellids
save(fullfile(resdir,fnmm),'phasic_cells','BetaIndex','GammaIndex')

% -------------------------------------------------------------------------
function [filenames]= sourcefiles(sourcedir)

files = dir(sourcedir);
filenames = cell(1,length(files));
for j = 3:length(files)
    filenames{j} = files(j).name;
end