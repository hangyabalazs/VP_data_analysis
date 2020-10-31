function [phasic_cells, BetaIndex, GammaIndex, isbeta, isgamma] = pr_freq_detection_VP(vpcells, resdir, resdir2, resdir3, acgsource, psthsource)
%PR_FREQ_DETECTION_VP finds rhythmic neurons based  on their autocorrelograms.
%
%   [PHASIC_CELLS, BETAINDEX, GAMMAINDEX, ISBETA, ISGAMMA] = PR_FREQ_DETECTION_VP(RESDIR,
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
phasewin1 = 10;
phasewin2 = 5;
ACGh = SCCR(:,1001:2000);  % half acg (smoothed)
lagsh = lags(1001:2000); % corresponding half of time vector in ms
for i = 1:length(cellids)  % loop through cells
    cACGh = ACGh(i,:);
    peaks(i) = max(cACGh);
    peak_inxx = find(cACGh == peaks(i));
    peak_inx(i)= peak_inxx(end);
    lagms(i) = lagsh(peak_inx(i));
end
phasic_cells = cellids(10<lagms & lagms<=67); % only choosing beta and gamma cells and eliminate bursting ones (lag>10ms)
phasic_lagms = lagms(10<lagms & lagms<=67);

% Beta and gamma index calculation
numPC = length(phasic_cells); % number of phasic cells
[betainx, gammainx, cellinx] = deal(nan(1,numPC));
[isgamma, isbeta] = deal(zeros(1,numPC));
for j = 1:numPC  % loop through non-bursting cells with <67 ms ACG peaks
    
    cellid = phasic_cells{j};
    cellinx(j) = find(strcmp(cellids,cellid)); % getting index from original cellidlist
    pACGh = ACGh(cellinx(j),:);  % acg of the current phasic cell
    
    % Find beta peak
    b_peak = pACGh(lagsh<=67 & lagsh>=34);
    [b_peak_val, binx] = max(b_peak);
    flank = phasewin1 / res;
    mbp = mean(pACGh(67+binx-flank:67+binx+flank)); % mean beta peak
    
    % Calculate baseline and Beta Index
    peak_loc = 67+binx;
    bsl = (pACGh(round(peak_loc/2))+pACGh(round(peak_loc*1.5))) / 2;  % baseline
    betainx(j) = (mbp-bsl)/max(mbp,bsl); % beta index
    
    % Find gamma peak
    g_peak = pACGh(lagsh<=33 & lagsh>=10);
    g_peak_inx = find(lagsh<=33 & lagsh>=10);
    [g_peak_val, ginx] = max(g_peak);
    flank = phasewin2 / res;
    mgp = mean(pACGh(max(19+ginx-flank,1):19+ginx+flank)); % mean gamma peak
    
    % Calculate baseline and Gamma Index
    peak_loc2 = 19 + ginx;
    bsl2 = (pACGh(round(peak_loc2/2))+pACGh(round(peak_loc2*1.5))) / 2;
    gammainx(j) = (mgp-bsl2)/max(mgp,bsl2); % gamma index
    
%     close all;figure;plot(lagsh,pACGh)
%     title([num2str(betainx(j)) ' ' num2str(gammainx(j))])
%     if betainx(j) > 0.4 || gammainx(j) > 0.25
%         hold on;plot(lagsh,pACGh,'r')
%     end
%     1;
end
BetaIndex = betainx;
GammaIndex = gammainx;
ph_Refractory = Refractory(cellinx);   % refractory period

% Gamma and beta cells
isbeta(setdiff(find(BetaIndex>0.4), find(GammaIndex>0.25))) = 1;
isgamma(setdiff(find(GammaIndex>0.25), find(BetaIndex>0.4))) = 1;

% Decide ambiguous cells
Ambig_cells = intersect(find(BetaIndex>0.4), find(GammaIndex>0.25));
for k = Ambig_cells
    if phasic_lagms(k) <= 67 && phasic_lagms(k) >= 34
        isbeta(k) = 1;
    elseif phasic_lagms(k) <= 33 && phasic_lagms(k) >= 10
        isgamma(k) = 1;
    end
end

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
save(fullfile(resdir,fnmm),'phasic_cells','BetaIndex','GammaIndex', 'isbeta', 'isgamma')

% -------------------------------------------------------------------------
function filenames = sourcefiles(sourcedir)

files = dir(sourcedir);
filenames = cell(1,length(files));
for j = 3:length(files)
    filenames{j} = files(j).name;
end