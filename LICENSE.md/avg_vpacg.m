function avg_vpacg(cellids,window,varargin)
%AVG_VPACG   Average autocorrelation.
%   AVG_VPACG loads ACGs and calculates average Z-scored ACG. For
%   Z-score-ing, the surrogate distribution is used (see
%   SOMCCG_CONF_FILTER).
%
%   See also VPACG and SOMCCG_CONF_FILTER.

%   Balazs Hangya, Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu
%   11-Nov-2018

%   Code review: BH 4/8/20

% Input arguments
prs = inputParser;
addRequired(prs,'cellids',@(s)iscell(s)|iscellstr(s)|ischar(s))
addRequired(prs,'window',@isscalar)  % ACG window, in seconds
addParameter(prs,'dt',0.0005,@isnumeric)   % ACG resolution, in seconds
addParameter(prs,'issave',false,@islogical)   % control saving behavior
addParameter(prs,'resdir',[],@(s)isdir(s)|isempty(s))   % results directory
addParameter(prs,'longsegments',false,@islogical)   % use only the longest segment after segment filtering
addParameter(prs,'seglim',0.3,@isnumeric);   % minimal segment length (s) to perform ACG calculation if 'longsegments' is 'true'
addParameter(prs,'segfilter','none',@(s)ischar(s)|iscellstr(s))   % filter segments
addParameter(prs,'filterinput',[])   % some filters need additional input
addParameter(prs,'minspikeno',100,@isnumeric)   % calculate ACG above minimal spike number
addParameter(prs,'maxspikeno',5000,@isnumeric)   % restrict included spikes
parse(prs,cellids,window,varargin{:})
g = prs.Results;
if ischar(cellids)
    cellids = {cellids};  % one cell ID
end

% Pass the control to the user in case of error
dbstop if error

% Determine time window
sr = 1000;      % sampling rate
wn = window * sr;    % CCG window in data points

% Cell loop for ACG
limit_spikes = [g.minspikeno g.maxspikeno];   % include max 50000 spikes; calculate only if min 100 spikes by default
numCells = length(cellids);
[CCR, LCCR, UCCR] = deal(zeros(numCells,2*wn));
[MeanH0, SDH0] = deal(nan(numCells,1));
SegmentLength = nan(numCells,1);

for iC = 1:numCells   % loop through the cells
    cell = cellids{iC};
    if isequal(g.segfilter,'none')
           ncc = loadcb(cell,'Spikes');   % use all spikes
    else
        tseg = findSegs3(cell,'segfilter',g.segfilter,...
            g.filterinput{:});  % find time segments
        ltseg = tseg(2,:) - tseg(1,:);  % length of the segments
        if g.longsegments   % use the longest segment if it's longer than the threshold ('seglim')
            seginx = find(ltseg==max(ltseg));
            tseg = tseg(:,seginx(1));   % find the longest segment
            ltseg = ltseg(seginx(1));
            if tseg < g.seglim * sr
                continue
            end
        end
        SegmentLength(iC) = sum(ltseg);  % cumulative length of the segments
        ncc = extractSegSpikes(cell,tseg);   % find spikes in the time segments
    end
    
    % Implement upper spike number limits
    if length(ncc) > limit_spikes(2)      % crop if too long to avoid out of memory
        ncc = ncc(1:limit_spikes(2));
    end
    
    if length(ncc) > limit_spikes(1)     % implement minimum number of spikes
        [H1, ccr, lwr, upr, rccg] = somccg_conf_filter(ncc,ncc,wn,1100);    % 1->2
        close(H1)

        CCR(iC,:) = ccr;   % auto-correlogram
        LCCR(iC,:) = lwr;  % lower significance limit
        UCCR(iC,:) = upr;  % upper significance limit
        MeanH0(iC) = mean(mean(rccg,2),1);   % surrogate mean
        SDH0(iC) = mean(std(rccg,[],2),1);   % surrogate SD
        disp(['Cell #' num2str(iC) ' / ' num2str(numCells) ' done......'])
    end
end

% Average z-scored ACG (by surrogate)
zCCR = nan(size(CCR));
for iP = 1:numCells
%     zCCR(iP,:) = (CCR(iP,:) - mean(CCR(iP,:))) / std(CCR(iP,:));
    zCCR(iP,:) = (CCR(iP,:) - MeanH0(iP)) / SDH0(iP);   % Z-score with surrogate mean and SD
end
avgccr = mean(zCCR);
serr = std(zCCR) / sqrt(numCells);

% Plot
H = figure;
time = linspace(-wn,wn,length(avgccr));
errorshade(time,avgccr,serr);
set(gca,'XLim',[-wn wn])