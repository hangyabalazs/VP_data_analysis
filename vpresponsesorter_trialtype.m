function vpresponsesorter_trialtype(cellids,issave,responsespec)
%VPRESPONSESORTER   Peri-event time histogram.
%   VPRESPONSESORTER(CELLIDS,ISSAVE) calculates non-adaptive PSTHs for a
%   set of cells (see ULTIMATE_PSTH) aligned to stimulus onset and feedback
%   delivery. Statistical tests are performed to probe significant firing
%   rate changes of responses (see PSTH_STATS). Indicators for significant
%   responses (p < 0.01) are added to CellBase as properties.
%   Input parameters:
%       CELLIDS - list of cell IDs; if empty or not specified, all
%           well-separated cells are selected (ID>20, L-ratio<0.15; see
%           LRATIO) from ventral pallidum
%       ISSAVE - controls saving
%       RESDIR - results directory
%       RESPONSESPEC - 'cue' for cue response, 'rew' for reward response
%       and 'pun' for punishment response
%
%   See also ULTIMATE_PSTH and LRATIO.

%   Balazs Hangya, Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu
%   11-Nov-2018

%   Code review: BH 7/24/19

% Directories
response_resdir = fullfile(resdir,'vpresponsesorter_trialtype');   % results directory
if ~isdir(response_resdir)
    mkdir(response_resdir)
end

% Input argument check
narginchk(0,4);
if nargin < 2
    issave = true;   % default saving behavior
end
if nargin < 1
    vpcells = [];   % all well-isolated units
else
    vpcells = vpselectcells(vpcells);
end
numCells = length(vpcells);

% Load CellBase
load(getpref('cellbase','fname'),'CELLIDLIST');

% Raster + PSTH
switch responsespec
    case 'cue'
        % Raster + PSTH aligned to stimulus onset
        alignevent = 'StimulusOn';   % trigger event
        shevent = 'DeliverAllFeedback';  % show-event
        partition = '#TrialType';   % partition trials
        wn = [-2 2];   % full raster window in seconds
        dt = 0.001;   % raster resolution, in seconds
        sigma = 0.02;   % controls smoothing for 'spsth'
        bwin = [-1 0];   % baseline window for MW-test
        twin = [0 0.5];   % test-window for MW-test
        

        % PSTH
        for iC = 1:numCells
            cellid = vpcells{iC};   % current cell
            stats1 = rasterPSTH(cellid,alignevent,shevent,partition,wn,dt,sigma,bwin,twin,issave,response_resdir);
        end
            if ~st
                error('MATLAB:vpresponsesorter:setvalueStatus','Error using setvalue.m')
            end
        end
        
    case 'rew'
        % Raster + PSTH aligned to feedback delivery
        alignevent = 'DeliverAllFeedback';   % trigger event
        shevent = 'StimulusOn';  % show-event
        partition = '#Reward==1';   % partition trials
        wn = [-3 3];   % full raster window in seconds
        dt = 0.001;   % raster resolution, in seconds
        sigma = 0.02;   % controls smoothing for 'spsth'
        bwin = [-3 -2];   % baseline window for MW-test (before cue) - -3.4 was the limit
        twin = [0 0.5];   % test-window for MW-test
        
        % PSTH
        for iC = 1:numCells
            cellid = vpcells{iC};   % current cell
            stats2 = rasterPSTH(cellid,alignevent,shevent,partition,wn,dt,sigma,bwin,twin,issave,response_resdir);
            end
            if ~st
                error('MATLAB:vpresponsesorter:setvalueStatus','Error using setvalue.m')
            end
            
            % Raster + PSTH aligned to feedback delivery
        alignevent = 'DeliverAllFeedback';   % trigger event
        shevent = 'StimulusOn';  % show-event
        partition = '#Reward==2';   % partition trials
        wn = [-3 3];   % full raster window in seconds
        dt = 0.001;   % raster resolution, in seconds
        sigma = 0.02;   % controls smoothing for 'spsth'
        bwin = [-3 -2];   % baseline window for MW-test (before cue) - -3.4 was the limit
        twin = [0 0.5];   % test-window for MW-test
        
        % PSTH
        for iC2 = 1:numCells
            cellid = vpcells{iC2};   % current cell
            stats3 = rasterPSTH(cellid,alignevent,shevent,partition,wn,dt,sigma,bwin,twin,issave,response_resdir);
            end
            if ~st
                error('MATLAB:vpresponsesorter:setvalueStatus','Error using setvalue.m')
            end
        end
        
        
    case 'pun'
        % Raster + PSTH aligned to feedback delivery
        alignevent = 'DeliverAllFeedback';   % trigger event
        shevent = 'StimulusOn';  % show-event
        partition = '#Punishment == 1';   % partition trials
        wn = [-3 3];   % full raster window in seconds
        dt = 0.001;   % raster resolution, in seconds
        sigma = 0.02;   % controls smoothing for 'spsth'
        bwin = [-3 -2];   % baseline window for MW-test (before cue) - -3.4 was the limit
        twin = [0 0.5];   % test-window for MW-test
        
        % PSTH
        for iC = 1:numCells
            cellid = vpcells{iC};   % current cell
            TE = loadcb(cellid,'TrialEvents');
            if ~all(isnan(TE.Punishment)) && isfield(TE, 'Punishment')  % if punishment was part of the session
                
                stats2 = rasterPSTH(cellid,alignevent,shevent,partition,wn,dt,sigma,bwin,twin,issave,response_resdir);
            end
        end
                
                if ~st
                    error('MATLAB:vpresponsesorter:setvalueStatus','Error using setvalue.m')
                end

% -------------------------------------------------------------------------
function stats1 = rasterPSTH(cellid,alignevent,shevent,partition,wn,dt,sigma,bwin,twin,issave,response_resdir)

% Raster plot and PSTH
TE = loadcb(cellid,'TrialEvents');   % load trial events
SP = loadcb(cellid,'EVENTSPIKES');
fld = fieldnames(TE);
if ~isequal(length(SP.event_stimes{1}),length(TE.(fld{1})))
    error('MATLAB:vpresponsesorter:rasterPSTH:trialMismatch',...
        'Trial number mismatch between TrialEvents and EVENTSPIKES.')
end

% Raster plot
viewcell2b(cellid,'TriggerName',alignevent,'SortEvent','TrialStart','sigma',sigma,...
    'eventtype','behav','ShowEvents',{{shevent}},'Partitions',partition,'window',wn,'PSTHPlot',false);
V_handle = gcf;
maximize_figure

% Peri-event time histogram
PSTHaxis_handle = findobj(allchild(V_handle),'type','axes','XLim',[0 1],'YLim',[0 1],'Tag','');   % handle for the empty PSTH axes
[~, ~, ~, ~, ~, stats1] = ...
    ultimate_psth(cellid,'trial',alignevent,wn,...
    'dt',dt,'sigma',sigma,'parts',partition,'isadaptive',0,...
    'maxtrialno',Inf,'baselinewin',bwin,'testwin',twin,'relative_threshold',0.1,'display',true); % calculate psth

% Plot & save
if ~iscell(stats1)
    stats1 = {stats1};   % only one PSTH
end
[~, tags] = partition_trials(TE,partition);
[mylabels, mycolors] = makeColorsLabels(@defineLabelsColors_Balazs,tags);
NumStats = length(stats1);
close_handle = nan(1,NumStats);
for iP =  1:NumStats
    Ls = findobj(allchild(stats1{iP}.axis_handle),'Type','line');   % all lines in the plot
    clr = findobj(Ls,'Type','line','Color','black');   % re-color the black PSTH
    set(clr,'Color',mycolors{iP})
    tx = findobj(allchild(stats1{iP}.axis_handle),'Type','text');   % re-position text object
    x_lim = xlim;
    y_lim = ylim;
    set(tx(1),'Position',[x_lim(1)+diff(x_lim)*0.1 (iP-1)*diff(y_lim)*0.4+y_lim(1)+diff(y_lim)*0.6 0])
    tx(1).String = regexprep(tx(1).String,'MW test',mylabels{iP});   % more meaningful labels
    if length(tx) > 1
        set(tx(2),'Position',[x_lim(1)+diff(x_lim)*0.1 (iP-1)*diff(y_lim)*0.4+y_lim(1)+diff(y_lim)*0.4  0])
        tx(2).String = regexprep(tx(2).String,'MW test',mylabels{iP});   % more meaningful labels
    end
    copyobj(tx,PSTHaxis_handle)   % copy PSTH and stat text to the raster figure
    copyobj(Ls,PSTHaxis_handle)
    xlabel(PSTHaxis_handle,['Time from ' alignevent])
    close_handle(iP) = stats1{iP}.figure_handle;
end

% Save figure
if issave
    cellidt = regexprep(cellid,'\.','_');
    fnm = fullfile(response_resdir,[cellidt '_' alignevent '_' partition(2:end) '.jpg']);
    saveas(V_handle,fnm)
    fnm2 = fullfile(response_resdir,[cellidt '_' alignevent '_' partition(2:end) '.fig']);
    saveas(V_handle,fnm2)
    fnm = fullfile(response_resdir,[cellidt '_' alignevent '_' partition(2:end) '.mat']);
    warning('off','MATLAB:Figure:FigureSavedToMATFile')
    save(fnm,'stats1')
end

close([V_handle, close_handle])