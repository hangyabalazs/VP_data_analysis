function lick_psth_summary_VP(vpcells,resdir,issave)
%LICK_PSTH_SUMMARY_VP   Average lick PETH.
%   LICK_PSTH_SUMMARY_VP plots average lick PETH (beam break time stamps
%   aligned to an event). 10 last sessions from mice are used.
%
%   See also ULTIMATE_PSTH.

%   Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hegedus.panna@koki.mta.hu
%   20-Dec-2019

%   Code review: BH 12/20/19, 1/2/20, 4/8/20

% Input arguments
if nargin < 1 || isempty(vpcells)
    vpcells = vpselectcells(vpcells);
end
if nargin < 2
    resdir = ('E:\auditory_pavlovian_cellbase\_paper_figs\code_revirew\Fig1_behavior\');   % result directory
end
if nargin < 3
    issave = true;
end

% Directories
if ~isfolder(resdir)
    mkdir(resdir)
end

% Select sessions: mice w no feedback delay that had enough training - last 10
% sessions (only first sessions if there were more sessions per day)
cells2use = [];
cellids = vpcells'; % VP cell IDs
animalIDs = getvalue('RatId', cellids); 
SessionID = getvalue('DateNum', cellids);
animals = unique(animalIDs); 
for i = 1:length(animals)  % find cells of all animals
    currentCells = cellids(animalIDs == animals(i));
    sessionIDs = sort(unique(getvalue('DateNum', currentCells)), 'ascend');
    if length(sessionIDs) >= 10
        sessions2use = sessionIDs(end-9:end);  % use only the last 10 sessions for the lickPSTH
    elseif length(sessionIDs) < 10
        sessions2use = sessionIDs;
    end
    for j =1:length(sessions2use)
        cellsofsession = cellids(animalIDs == animals(i) & SessionID == sessions2use(j));
        [aID, sID] = cellid2tags(cellsofsession{1});
        cells2use = [cells2use; {aID sID}];   %#ok<*CCAT1> % take one cell per session
    end
end

% Time window
wn = [-5 5];   % in seconds
dt = 0.001;   % resolution, in seconds
time = wn(1):dt:wn(2);   % time vector

% PETH
NumSessions = size(cells2use',2);   % number of sessions
[Hit_allpsth, FA_allpsth] = deal([]);
for iS = 1:NumSessions
    sessionid = cells2use(iS,:);
    [spsth_hit, spsth_fa] = main(sessionid,wn,dt);
    Hit_allpsth = [Hit_allpsth; spsth_hit]; %#ok<*AGROW>
    FA_allpsth = [FA_allpsth; spsth_fa]; %#ok<*AGROW>
    H = gcf;
    if issave
        cellidt = [sessionid{1} '_' sessionid{2}];
        fnm = fullfile(resdir, [cellidt '_LICK.fig']);   % save
        saveas(H,fnm)
        fnm = fullfile(resdir, [cellidt '_LICK.jpg']);
        saveas(H,fnm)
    end
    close(H)
end

% Plot & save
H = figure;
green = [51 204 51] / 255;   % colors for plotting
red = [216 41 0] / 255;
errorshade(time,nanmean(Hit_allpsth),nanstd(Hit_allpsth)/sqrt(size(Hit_allpsth,1)),...
    'LineColor',green,'ShadeColor',green)
hold on
errorshade(time,nanmean(FA_allpsth),nanstd(FA_allpsth)/sqrt(size(FA_allpsth,1)),...
    'LineColor',red,'ShadeColor',red)
if issave
    fnm = fullfile(resdir,'average_lickPSTH.fig');
    saveas(H,fnm)
end

% -------------------------------------------------------------------------
function [spsth_hit, spsth_fa] = main(cellid,wn,dt)

% Filter input
filterinput_hit = 'TrialType==1';
filterinput_fa = 'TrialType==2';

% Calcualte lick PSTH
[~, spsth_hit] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.02,'event_filter','custom','filterinput',filterinput_hit,...
    'isadaptive',0,'maxtrialno',Inf);
[~, spsth_fa] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.02,'event_filter','custom','filterinput',filterinput_fa,...
    'isadaptive',0,'maxtrialno',Inf);

% Lick raster
H = figure;
viewlick(cellid,'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav',...
    'ShowEvents',{{'StimulusOn' 'StimulusOff'}},...
    'Partitions','#TrialType','window',wn)
maximize_figure(H)