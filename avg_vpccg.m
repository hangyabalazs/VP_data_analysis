function avg_vpccg
%AVG_VPCCG   Average cross-correlation.
%
%   AVG_VPCCG loads CCG pairs and calculates average Z-scored CCG. For
%   Z-score-ing, the surrogate distribution is used (see
%   SOMCCG_CONF_FILTER).
%
%   See also VPCCG and SOMCCG_CONF_FILTER.

%   Balazs Hangya, Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu
%   11-Nov-2018

%   Code review: BH 2/12/20

% Load pairs
load([getpref('cellbase','datapath') '\VP\vpccg\cellgroups_tetrodepairs.mat']);
load([getpref('cellbase','datapath') '\VP\vpccg\cellgroups_nontetrodepairs.mat']);
numPairs = length(sync_exc_nttp)+length(sync_exc_ttp)+length(sync_monosyn_nttp); %number of Pairs
PairOfCells = [sync_exc_nttp sync_exc_ttp sync_monosyn_nttp];

% Parameters
segfilter = 'stim_excl_vp';
filterinput = {'light_activation_duration',[-5 5],'margins',[0 0]};
window = 0.05;
minspikeno = 100;
maxspikeno = 10000;
longsegments = false;
seglim = 0.3;

% Determine time window
sr = 1000;      % sampling rate
wn = window * sr;    % CCG window in data points

% Cell pair loop for CCG
wb = waitbar(0,'Please wait...','Name','Running CCG...');  % progress indicator
global WB
WB(end+1) = wb;
limit_spikes = [minspikeno maxspikeno];   % include max 50000 spikes; calculate only if min 100 spikes by default
[CCR, LCCR, UCCR] = deal(zeros(numPairs,2*wn+1));
[MeanH0, SDH0] = deal(nan(numPairs,1));
SegmentLength = nan(numPairs,1);
for iP = 1:numPairs   % loop through pairs of cells
    cell1 = PairOfCells{iP}{1};
    cell2 = PairOfCells{iP}{2};
    if isequal(segfilter,'none')
        ncc1 = loadcb(cell1,'SPIKES');   % use all spikes
        ncc2 = loadcb(cell2,'SPIKES');
    else
        tseg = findSegs3(cell1,'segfilter',segfilter,...
            filterinput{:});  % find time segments
        ltseg = tseg(2,:) - tseg(1,:);  % length of the segments
        if longsegments   % use the longest segment if it's longer than the threshold ('seglim')
            seginx = find(ltseg==max(ltseg));
            tseg = tseg(:,seginx(1));   % find the longest segment
            ltseg = ltseg(seginx(1));
            if tseg < seglim * sr
                continue
            end
        end
        SegmentLength(iP) = sum(ltseg);  % cumulative length of the segments
        ncc1 = extractSegSpikes(cell1,tseg);   % find spikes in the time segments
        ncc2 = extractSegSpikes(cell2,tseg);
    end
    
    % Implement upper spike number limits
    if length(ncc1) > limit_spikes(2)      % crop if too long to avoid out of memory
        ncc1 = ncc1(1:limit_spikes(2));
    end
    if length(ncc2) > limit_spikes(2)
        ncc2 = ncc2(1:limit_spikes(2));
    end
    
    if length(ncc1) > limit_spikes(1) && length(ncc2) > limit_spikes(1)     % implement minimum number of spikes
        [H1, ccr, lwr, upr, rccg] = somccg_conf_filter(ncc1,ncc2,wn,1100);    % 1->2
        xl = xlim;   % put the cell IDs in the plot
        yl = ylim;
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.9,regexprep(cell1,'_',' '))
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.85,regexprep(cell2,'_',' '))
        close(H1)
        
        CCR(iP,:) = ccr;   % cross-correlogram
        LCCR(iP,:) = lwr;  % lower significance limit
        UCCR(iP,:) = upr;  % upper significance limit
        MeanH0(iP) = mean(mean(rccg,2),1);   % surrogate mean
        SDH0(iP) = mean(std(rccg,[],2),1);   % surrogate SD
    end
end

% Average z-scored CCG (by surrogate)
zCCR = nan(size(CCR));
for iP = 1:numPairs
%     zCCR(iP,:) = (CCR(iP,:) - mean(CCR(iP,:))) / std(CCR(iP,:));
    zCCR(iP,:) = (CCR(iP,:) - MeanH0(iP)) / SDH0(iP);   % Z-score with surrogate mean and SD
end
avgccr = mean(zCCR);
serr = std(zCCR) / sqrt(numPairs);

% Plot
H1 = figure;
time = linspace(-wn,wn,length(avgccr));
errorshade(time,avgccr,serr);
set(gca,'XLim',[-wn wn])