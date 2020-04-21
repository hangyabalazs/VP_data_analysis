function compare_expectations_VP(vpcells, responsetype, responsespec, resdir)
%COMPARE_EXPECTATIONS_VP   Average, normalized peri-event time-histograms
%comparing neuronal response to expected and unexpected stimuli.
%
%   COMPARE_EXPECTATIONS_VP(VPCELLS, RESPONSESPEC, RESDIR) Calculates
%   normalized average PSTHs  of cells defined by VPCELLS, aligned to
%   RESPONSETYPE (cue, reward or punishment) and compares neuronal response
%   to expected and unexpected stimuli (timestamps extracted from
%   TrialEvents). Both activation (RESPONSESPEC = 1) and inhibition
%   (RESPONSESPEC = -1) to the stimulus are compared. The PSTHs are saved
%   to RESDIR.
%
%   See also ULTIMATE_PSTH.

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   05-Feb-2020 

%   Code review: BH 2/10/20, 4/8/20

% Input arguments
if nargin < 1 || isempty(vpcells)
    vpcells = vpselectcells(vpcells);
end

% Directories
if ~isfolder(resdir)
    mkdir(resdir)
end
    
% Sub-select neurons that respond to the event of interest
resp = getvalue(responsetype,vpcells);
cells = vpcells(resp==responsespec);

% Exclude those cells, where there is only one trialtype (or low number of trials)
NumCells = length(cells);
delstr = ones(1,NumCells);
for i = 1:NumCells
    cellid = cells(i);
    E = loadcb(cellid,'TrialEvents');   % load events
    if strcmp(responsetype, 'cueresponse')
        if sum(E.TrialType==2) == 0  % if there is no second type of trials
             delstr(i) = 0;
        end
    elseif strcmp(responsetype, 'rewardresponse')
        if (sum(E.AllReward==2) == 0) || (sum(E.AllReward==2) < 5) || (sum(E.AllReward==1) < 5) % if there is no second type of trials
            delstr(i) = 0;
        end
    elseif strcmp(responsetype, 'punishresponse')
        if (sum(E.Punishment==1) == 0) || (sum(E.Punishment==2) < 5) || (sum(E.Punishment==1) < 5) % if there is no second type of trials
            delstr(i) = 0;
        end
    end
end
cells = cells(delstr==1);

% Define filterinput
switch responsetype
    case 'cueresponse'
        filters = {'TrialType==1' 'TrialType==2'};
        align = 'StimulusOn';
    case 'rewardresponse'
        filters = {'AllReward==1' 'AllReward==2'};
        align = 'DeliverAllFeedback';
    case 'punishresponse'
        filters = {'Punishment==1' 'Punishment==2'};
        align = 'DeliverAllFeedback';
end

% Time window
wn = [-4 4];
dt = 0.001;
time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector

% PSTH
R2 = runanalysis(@ultimate_psth,...
    'trial', align, wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
    'event_filter','custom', 'filterinput',filters{2},'maxtrialno',Inf,...
    'baselinewin',[-3 -2],'testwin',[0 0.5],'relative_threshold',0.1,'cellids',cells);
R1 = runanalysis(@ultimate_psth,...
    'trial', align, wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
    'event_filter','custom', 'filterinput',filters{1},'maxtrialno',Inf,...
    'baselinewin',[-3 -2],'testwin',[0 0.5],'relative_threshold',0.1,'cellids',cells);

% Z-score by baseline
NumCell = length(cells);
[psth_R1, psth_R2] = deal(nan(NumCell,size(R1{1,1},2)));
baseline_inx = 1:800;
Rinx = 2;
for iC = 1:NumCell
    mn = mean(R1{iC,Rinx}(baseline_inx));
    sd = std(R1{iC,Rinx}(baseline_inx));
    sd2 = std(R2{iC,Rinx}(baseline_inx));
    
    if sd >= 0.1 && sd2 >= 0.1
        psth_R1(iC,:) = (R1{iC,Rinx} - mn) / sd;   % the same mean and SD used within a cell,
        psth_R2(iC,:) = (R2{iC,Rinx} - mn) / sd;   % to preserve the trial-type differences
    end
end
psth_R1(any(isnan(psth_R1), 2), :) = [];
psth_R2(any(isnan(psth_R2), 2), :) = [];

% Plot
inx = 1:size((psth_R1),1);  % all
disp(length(inx))
H = figure;
errorshade(time,mean(psth_R1(inx,:)),std(psth_R1(inx,:))/sqrt(size(psth_R1(inx,:),1)),...
    'LineColor',[0 0.9 0],'ShadeColor',[0 0.9 0])
errorshade(time,mean(psth_R2(inx,:)),std(psth_R2(inx,:))/sqrt(size(psth_R2(inx,:),1)),...
    'LineColor',[0.5 0 0],'ShadeColor',[0.5 0 0])
set(gca,'XLim',[-3000 3000]);

% Save figure
if responsespec ==1
    tag = 'excitation';
elseif responsespec == -1
    tag = 'inhibition';
end
filename = fullfile(resdir,['compare_expectation' responsetype '_' tag '.fig']);
saveas(H,filename);