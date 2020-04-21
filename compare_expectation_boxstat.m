function compare_expectation_boxstat(vpcells, responsetype, responsespec, resdir)
%COMPARE_EXPECTATION_BOXSTAT compares neuronal response to expected and
%unexpected stimuli with a Wilcoxon signed rank test.
%
%   COMPARE_EXPECTATION_BOXSTAT(VPCELLS, RESPONSETYPE, RESPONSESPEC,
%   RESDIR) compares neuronal response to expected and unexpected stimuli
%   (RESPONSETYPE) with a Wilcoxon signed rank test. Cellids are listed in
%   VPCELLS. Both activation (RESPONSESPEC =1) and inhibition (RESPONSESPEC
%   = -1) to the reward and punishment predicting stimului are compared.
%   The results are saved to RESDIR as box-whiskers plots.
%
%   See also ULTIMATE_PSTH and BOXSTAT.

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
        if sum(E.TrialType==2) == 0 
             delstr(i) = 0;
        end
    elseif strcmp(responsetype, 'rewardresponse')
        if (sum(E.AllReward==2) == 0) || (sum(E.AllReward==2) < 5) || (sum(E.AllReward==1) < 5)  % if there is no second type of trials
            delstr(i) = 0;
        end
    elseif strcmp(responsetype, 'punishresponse')
        if (sum(E.Punishment==1) == 0) || (sum(E.Punishment==2) < 5) || (sum(E.Punishment==1) < 5)  % if there is no second type of trials
            delstr(i) = 0;
        end
    end
end
cells = cells(delstr==1);

% Time window
wn = [-4 4];
dt = 0.001;
time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector

% PSTH
R2 = runanalysis(@ultimate_psth,...
    'trial', 'StimulusOn', wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
    'event_filter','custom', 'filterinput','TrialType==2','maxtrialno',Inf,...
    'baselinewin',[-1 0],'testwin',[0 0.5],'relative_threshold',0.1,'cellids',cells);
R1 = runanalysis(@ultimate_psth,...
    'trial', 'StimulusOn', wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
    'event_filter','custom', 'filterinput','TrialType==1','maxtrialno',Inf,...
    'baselinewin',[-1 0],'testwin',[0 0.5],'relative_threshold',0.1,'cellids',cells);

% Index for time 0
bwin = [-3 -2];   % baseline window for MW-test
twin = [0 0.5];   % test-window for MW-test
st = abs(wn(1)) / dt;   % in ms
nullindex = st + 1;

% Window for testing the potential effect
WNb = [bwin(1)/dt+nullindex bwin(2)/dt+nullindex-1];   % baseline window; convert to indices; exlude 0
WNt = [twin(1)/dt+nullindex twin(2)/dt+nullindex-1];   % test window; convert to indices
WNb = round(WNb);
WNt = round(WNt);
lWNb = WNb(2) - WNb(1) + 1;   % length of baseline window
lWNt = WNt(2) - WNt(1) + 1;   % length of test window

BR1 = R1(:,5); % binraster for cue1
BR2 = R2(:,5); % binraster for cue2

% Spike counts (test window - baseline window)
avg_spikecount1 = nan(1, size(BR1 ,1));
avg_spikecount2 = nan(1, size(BR2 ,1));
for i =1:size(BR1 ,1)  % loop through cells
    avg_spikecount1_cell = nan(1,size(BR1{i} ,1));
    avg_spikecount2_cell = nan(1,size(BR2{i} ,1));
    for j = 1:size(BR1{i} ,1)   % loop through trials
        avg_spikecount1_cell(j) = sum(BR1{i}(j, WNt(1):WNt(2)), 2) / lWNt / dt - ...
            sum(BR1{i}(j, WNb(1):WNb(2))) / lWNb / dt;
    end
    avg_spikecount1(i) = mean(avg_spikecount1_cell);
    
    for k = 1:size(BR2{i} ,1)   % loop through trials
        avg_spikecount2_cell(k) = sum(BR2{i}(k, WNt(1):WNt(2)), 2) / lWNt / dt - ...
            sum(BR2{i}(k, WNb(1):WNb(2))) / lWNb / dt;
    end
    avg_spikecount2(i) = mean(avg_spikecount2_cell);
end

% Box-whisker plot
[H, Wp] = boxstat(avg_spikecount1,avg_spikecount2,'T1','T2',0.005,'paired');

% Save figure
if responsespec == 1
    tag = 'excitation';
elseif responsespec == -1
    tag = 'inhibition';
end
filename = fullfile(resdir,['compare_expectation_' responsetype '_' tag '_boxstat.fig']);
saveas(H,filename);