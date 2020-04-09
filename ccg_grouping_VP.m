function ccg_grouping_VP(cellids,varargin)
%CCG_GROUPING_VP  group cell pairs based on their ccg.
%   CCG_GROUPING_VP calculates cross-correlations. Window size is set to +-50 ms with
%   a 1 ms resolution. Maximum 50000 spikes are included to avoid memory
%   problems. CCG is not calculated if one of the cells has less than 100
%   spikes. Segments are not filtered. Minimal shift for
%   shuffled CCGs is set to 1100 ms. For details on the algorithm, see
%   SOMCCG_CONF_FILTER.
%
%   CCG_GROUPING_VP(I) calls CCG for pairs of cells defined by the cell ID list
%   (or index set to CELLIDLIST, see CellBase documentation) I. By default,
%   non-tetrode pairs within the given list are selected for analysis.
%   Optional input arguments (parameter-value pairs with default values):
%       'issave', false - controls saving behavior; plots and
%           cross-correlation matrices with confidence intervals are saved
%           only if 'issave' is set to true
%       'whichcells', 'nontetrodepairs' - method of pair selection;
%           'nontetrodepairs' selects cells from other tetrodes,
%           'tetrodepairs' selects cells from the same tetrode and
%           'allpairs' selects all cells from the session
%       'include', 'list' - by default, only pairs for which both cells are
%           included in I are analyzed; if 'include' is set to 'cellbase',
%           all cells in CellBase that are paired with the ones in I
%           according to 'whichcells' are analyzed
%
%   See also ACG, SOMCCG_CONF_FILTER and XCORR_WRAND_FILTER.

%   Balazs Hangya, Panna Hegedus 2018
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hegedus.panna@koki.mta.hu

% Input arguments
prs = inputParser;
addRequired(prs,'cellids',@(s)iscell(s)|iscellstr(s)|ischar(s))
% addRequired(prs,'window',@isscalar)  % CCG window, in seconds
addParameter(prs,'issave',false,@islogical)   % control saving behavior
addParameter(prs,'resdir',[],@(s)isdir(s)|isempty(s))   % results directory
addParameter(prs,'longsegments',false,@islogical)   % use only the longest segment after segment filtering
addParameter(prs,'seglim',0.3,@isnumeric);   % minimal segment length (s) to perform CCG calculation if 'longsegments' is 'true'
addParameter(prs,'segfilter','none',@(s)ischar(s)|iscellstr(s))   % filter segments
addParameter(prs,'filterinput',[])   % some filters need additional input
addParameter(prs,'minspikeno',100,@isnumeric)   % calculate CCG above minimal spike number
addParameter(prs,'maxspikeno',5000,@isnumeric)   % maximize included spikes
addParameter(prs,'whichcells','nontetrodepairs',...
    @(s)ismember(s,{'nontetrodepairs','tetrodepairs','allpairs'}))   % which cells to include
addParameter(prs,'include','list',...
    @(s)ismember(s,{'list','cellbase'}))   % cell pair selection behavior: only from input list or full CellBase
parse(prs,cellids,varargin{:})
g = prs.Results;

% Pass the control to the user in case of error
dbstop if error

% Directories
resdir1 = [g.resdir '\sorted_ccg' ]; % directory for acgs and ccgs
resdir2 = [resdir1 '\sync_monosyn_exc' ]; % directory for acgs and ccgs
resdir3 = [resdir1 '\sync_exc' ]; % directory for acgs and ccgs
resdir4 = [resdir1 '\monosyn_exc' ]; % directory for acgs and ccgs

if ~isfolder(resdir1)
    mkdir(resdir1);
end
if ~isfolder(resdir2)
    mkdir(resdir2);
end
if ~isfolder(resdir3)
    mkdir(resdir3);
end
if ~isfolder(resdir4)
    mkdir(resdir4);
end

% Determine time window
sr = 1000;      % sampling rate
wn = 50 * sr / 1000;    % 2 * 1000 ms window

% Cell pairs
PairOfCells = cell(0,2);
numCells = length(cellids);  % number of cells
for iC = 1:numCells
    cellid = cellids{iC};
    [animalID, sessionID, tetrode1, unit1] = cellid2tags(cellid);
    switch g.whichcells
        case 'nontetrodepairs'
            [nm, ps] = nontetrodepairs(cellid);   % cells on other tetrodes
        case 'tetrodepairs'
            [~, ps] = tetrodepairs(cellid);   % cells on the same tetrode
            ps = setdiff(ps,cellid);
            nm = length(ps);
        case 'allpairs'
            ps = findcell('rat',animalID,'session',sessionID);   % all concurrently recorded cells
            ps = setdiff(ps,cellid);
            nm = length(ps);
    end
    if isequal(g.include,'list')
        ps = intersect(cellids,ps);   % include only those that are in the input list
        nm = length(ps);
    end
    for k = 1:nm
        [~, ~, tetrode2, unit2] = cellid2tags(ps(k));
        if (tetrode2*10+unit2) > (tetrode1*10+unit1)   % prevent duplicates
            PairOfCells(end+1,1:2) = {cellid ps{k}}; %#ok<AGROW>
        end
    end
end
numPairs = size(PairOfCells,1);

% Cell pair loop for CCG
wb = waitbar(0,'Please wait...','Name','Running CCG...');  % progress indicator
global WB
WB(end+1) = wb;
limit_spikes = [g.minspikeno g.maxspikeno];   % include max 50000 spikes; calculate only if min 100 spikes
[CCR, LCCR, UCCR] = deal(zeros(numPairs,2*wn+1));
[MeanH0, SDH0] = deal(nan(numPairs,1));
SegmentLength = nan(numPairs,1);

[inx1, inx2, inx3, inx4, inx5] = deal(1);
[monosyn_exc_ttp, sync_exc_ttp, monosyn_exc_nttp, sync_exc_nttp, sync_monosyn_nttp] = deal({});
for iP = 1:numPairs   % loop through pairs of cells
    cell1 = PairOfCells{iP,1};
    cell2 = PairOfCells{iP,2};
    
    ncl1 = regexprep(cell1,'\.','_');
    ncl2 = regexprep(cell2,'\.','_');
    
    if isequal(g.segfilter,'none')
        ncc1 = loadcb(cell1,'SPIKES');   % use all spikes
        ncc2 = loadcb(cell2,'SPIKES');
    else
        tseg = findSegs3(cell1,'segfilter',g.segfilter,...
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
    
    if length(ncc1) > limit_spikes(1) && length(ncc2) > limit_spikes(1)     % minimum 100 spikes
        [H1 ccr lwr upr rccg] = somccg_conf_filter(ncc1,ncc2,wn,1100);    % 1->2
        %         close(H1)
    end
    
    maxval = max(ccr);
    maxinx = find(ccr == maxval);
    
    if (maxval-upr(maxinx)) > 10 % a better criterium for excitation may be defined later
        einx = find(ccr>upr); % find the indices, where there were excitation
        if isempty(einx)
            einx = 0;
        end
        
        if length(einx) > 1 % how 'wide' is the excitation around the peak
            segs = diff(einx);
            width = sum(segs==1) + 1;
        elseif length(einx) == 1
            width = 1;
        else
            width = 0;
        end
        fnm2 = ['CCG_' ncl1 '_' ncl2 '.fig'];
        fnm3 = ['CCG_' ncl1 '_' ncl2 '.jpg'];
        switch g.whichcells
            case 'nontetrodepairs' % Examining nontetrodepairs
                if width < 3 % 'Narrow' excitation
                    if max(einx) > 46 && max(einx) < 56
                        if ~isnan(find(einx==51))
                            sync_monosyn_nttp{inx1} = [{cell1}, {cell2}];   % narrow sync activation
                            inx1 = inx1 + 1;
                            saveas(H1,fullfile(resdir2,fnm2));   % save CCG plot
                            saveas(H1,fullfile(resdir2,fnm3));   % save CCG plot
                        else
                            monosyn_exc_nttp{inx2} = [{cell1}, {cell2}];    % monosynaptic excitation
                            inx2 = inx2 + 1;
                            saveas(H1,fullfile(resdir4,fnm2));   % save CCG plot
                            saveas(H1,fullfile(resdir4,fnm3));   % save CCG plot
                        end
                    end
                    
                    
                elseif width >= 3 % 'Broader' excitation
                    sync_exc_nttp{inx3} = [{cell1}, cell2];   % broad sync excitation
                    inx3 = inx3 + 1;
                    saveas(H1,fullfile(resdir3,fnm2));   % save CCG plot
                    saveas(H1,fullfile(resdir3,fnm3));   % save CCG plot
                end
                
            case 'tetrodepairs' % Examining tetrodepairs
                if width < 3
                    if max(einx) > 46 && max(einx) < 56
                        if ~isnan(find(einx==51)) % if there is a peak in the 0 bin = shared spikes, ignore these pairs
                            continue
                        else % if the peak is not right in the middle
                            monosyn_exc_ttp{inx4} = [{cell1}, {cell2}];  % monosynaptic exxcitation
                            inx4 = inx4 + 1;
                            saveas(H1,fullfile(resdir4,fnm2));   % save CCG plot
                            saveas(H1,fullfile(resdir4,fnm3));   % save CCG plot
                        end
                    end 
                elseif width >= 3
                    sync_exc_ttp{inx5} = [{cell1}, cell2];   % broad sync excitation
                    inx5 = inx5 + 1;
                    saveas(H1,fullfile(resdir3,fnm2));   % save CCG plot
                    saveas(H1,fullfile(resdir3,fnm3));   % save CCG plot
                    close(H1)
                end
        end
    else
        saveas(H1,fullfile(resdir1,fnm2));   % save CCG plot
        saveas(H1,fullfile(resdir1,fnm3));   % save CCG plot
    end
    close all
end
switch g.whichcells
    case 'nontetrodepairs'
        save(fullfile(resdir1,'cellgroups_nontetrodepairs.mat'),'monosyn_exc_nttp','sync_exc_nttp','sync_monosyn_nttp');
    case 'tetrodepairs'
        save(fullfile(resdir1,'cellgroups_tetrodepairs.mat'),'monosyn_exc_ttp','sync_exc_ttp');
end
end
