function extract_baselineFR(cellids,responsespec, resdir)
%EXTRACT_BASELINEFR   Calculates baseline firing rate of neurons.
%   EXTRACT_BASELINEFR(CELLIDS, RESPONSESPEC, RESDIR) calculates baseline
%   firing rate of neurons (CELLIDS) from a binraster aligned to
%   RESPONSESPEC. Results are saved to RESDIR and added to CellBase.
%
%   See also ULTIMATE_PSTH

%   Balazs Hangya, Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu
%   18-02-2021


if ~isfolder(resdir) % make results directory
    mkdir(resdir);
end

dbstop if error

% Aligning: cue, reward or punishment
switch responsespec
    case {'cue'}
        alignevent = 'StimulusOn';
        filterevent = 'StimulusOn';
        bslnwn = [-1 0]; % baseline window
    case {'reward'}
        alignevent = 'DeliverAllFeedback';
        filterevent = 'RewardedTrials==1';
        bslnwn = [-3 -2]; % baseline window
    case {'punishment'}
        alignevent = 'DeliverAllFeedback';
        filterevent = 'PunishedTrials==1';
        bslnwn = [-3 -2]; % baseline window
    otherwise
        error('MATLAB:avg_psth_VP:responsespecUnknown','Unrecognized input value for responsespec.')
end

% Load CELLIDLIST
load(getpref('cellbase','fname'),'CELLIDLIST');

if ~ismember('baseline_FR',listtag('prop')) % Add baseline firing rate property
    insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name','baseline_FR')
end

% Time windows
wn = [-4 4];  % full data window
dt = 0.001;  % time resolution
bsl_length = diff(bslnwn);   % baseline window length in s
time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector in ms
baseline_inx = find(time == bslnwn(1)*1000):find(time == bslnwn(2)*1000);   % baseline window indices wrt. to full data window

% Calculate firing rates
    NumCells = length(cellids);   % number of analyzed cells
    BasFR = nan(1,length(NumCells)); % preallocate space for baseline fr
    
    for i = 1:NumCells   % loop through cells
        cellid = cellids{i}; % current cell
        
        % Calculate bin raster with 'ultimate_psth'
        [~, ~, ~, ~, spt, ~] = ultimate_psth( cellid, 'trial', alignevent, wn,...
            'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
            'event_filter','custom', 'filterinput', filterevent,'maxtrialno',Inf,...
            'relative_threshold',0.1);
        
        numTrials = size(spt,1);  % number of trials
        baseline = nan(1,numTrials);
        for k = 1:numTrials
            baseline(k) = sum(spt(k, baseline_inx)) / bsl_length;   % firing rate in baseline window
        end
        BasFR(1,i) = mean(baseline);  % store firing rates
        st = setvalue(cellid,'baseline_FR',BasFR(1,i));   % add data to CellBase
    end
    
save(fullfile(resdir, ['compare_baseline_' responsespec '.mat']), 'BasFR') % save data