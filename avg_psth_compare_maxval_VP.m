function avg_psth_compare_maxval_VP(cellids1, responsespec1, cellids2, responsespec2, tag, resdir)
%AVG_PSTH_COMPARE_MAXVAL_VP   Comparison of reinforcement response.
%   AVG_PSTH_COMPARE_MAXVAL_VP(CELLIDS1,RESPONSESPEC1, CELLIDS2, RESPONSESPEC2, RESDIR) calculates 
%   average PSTHs and compares neuronal response to reward and punishment
%   (RESPONSESPEC) with a Mann-Whitney U-test. Both activation and inhibition 
%   can be compared (indicated by TAG). Cellids are listed in
%   CELLIDS1 and CELLIDS2 for reward and punishment response, respectively.
%   The results are saved to RESDIR as box-whiskers plots.
%
%   See also ULTIMATE_PSTH.

%   Panna Hegedus
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   03-Dec-2019

%   Code review: BH 4/8/20

% Directories
if ~isfolder(resdir)
    mkdir(resdir)
end

% Aligning: cue, reward or punishment
switch responsespec1
    case {'cue','cueresponse'}
        alignevent1 = 'StimulusOn';
        filterevent1 = 'StimulusOn';
    case {'reward','rewardresponse'}
        alignevent1 = 'DeliverAllFeedback';
        filterevent1 = 'RewardedTrials==1';
    case {'punishment','punishresponse','punishmentresponse'}
        alignevent1 = 'DeliverAllFeedback';
        filterevent1 = 'PunishedTrials==1';
    otherwise
        error('MATLAB:avg_psth_compare_maxval:responsespecUnknown','Unrecognized input value for responsespec.')
end

switch responsespec2
    case {'cue','cueresponse'}
        alignevent2 = 'StimulusOn';
        filterevent2 = 'StimulusOn';
        bwin = [-1 0]; % baseline window
    case {'reward','rewardresponse'}
        alignevent2 = 'DeliverAllFeedback';
        filterevent2 = 'RewardedTrials==1';
        bwin = [-3 -2]; % baseline window
    case {'punishment','punishresponse','punishmentresponse'}
        alignevent2 = 'DeliverAllFeedback';
        filterevent2 = 'PunishedTrials==1';
        bwin = [-3 -2]; % baseline window
    otherwise
        error('MATLAB:avg_psth_compare_maxval:responsespecUnknown','Unrecognized input value for responsespec.')
end

% Time
wn = [-4 4]; % time window
dt = 0.001; % time resolution
twin = [0 0.5];   % test-window for MW-test
st = abs(wn(1)) / dt;
nullindex = st + 1;   % index for time 0
% time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector
% btime = bwin(1)*1000:dt*1000:bwin(2)*1000;   % baseline time vector
% [~, baseline_inx] = intersect(time,btime);   % baseline indices

% Ultimate PSTH
R1 = runanalysis(@ultimate_psth,...
    'trial', alignevent1, wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
    'event_filter','custom', 'filterinput',filterevent1,'maxtrialno',Inf,...
    'baselinewin',bwin,'testwin',twin,'relative_threshold',0.1,'cellids',cellids1);
R2 = runanalysis(@ultimate_psth,...
    'trial', alignevent2, wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
    'event_filter','custom', 'filterinput',filterevent2,'maxtrialno',Inf,...
    'baselinewin',bwin,'testwin',twin,'relative_threshold',0.1,'cellids',cellids2);

% Window for testing the potential effect
WNb = [bwin(1)/dt+nullindex bwin(2)/dt+nullindex-1];   % baseline window; convert to indices; exclude 0
WNt = [twin(1)/dt+nullindex twin(2)/dt+nullindex-1];   % test window; convert to indices
WNb = round(WNb);
WNt = round(WNt);
lWNb = WNb(2) - WNb(1) + 1;   % length of baseline window
lWNt = WNt(2) - WNt(1) + 1;   % length of test window

BR1 = R1(:,5); % binraster 1
BR2 = R2(:,5); % binraster 2

% Firing rates (test window - baseline window)
avg_spikecount1 = nan(1, size(BR1 ,1));
avg_spikecount2 = nan(1, size(BR2 ,1));
for i = 1:size(BR1 ,1)  % loop through cells
    avg_spikecount1_cell = nan(1,size(BR1{i} ,1));
    for j = 1:size(BR1{i} ,1)   % loop through trials
        avg_spikecount1_cell(j) = sum(BR1{i}(j, WNt(1):WNt(2)), 2) / lWNt / dt - ...
            sum(BR1{i}(j, WNb(1):WNb(2))) / lWNb / dt;
    end
    avg_spikecount1(i) = mean(avg_spikecount1_cell);
end

for i = 1:size(BR2 ,1)  % loop through cells
        avg_spikecount2_cell = nan(1,size(BR2{i} ,1));
    for k = 1:size(BR2{i} ,1)   % loop through trials
        avg_spikecount2_cell(k) = sum(BR2{i}(k, WNt(1):WNt(2)), 2) / lWNt / dt - ...
            sum(BR2{i}(k, WNb(1):WNb(2))) / lWNb / dt;
    end
    avg_spikecount2(i) = mean(avg_spikecount2_cell);
end

% Box-whisker plot
[H, Wp] = boxstat(avg_spikecount1, avg_spikecount2,'T1','T2', 0.005); %#ok<ASGLU>
h = findobj(gca,'Type','line');
set(h,'Marker','none');

% Save figure
filename = fullfile(resdir,['compare_maxvalue_' responsespec1 '_' responsespec2 '_' tag  '_boxstat.fig']);
filename2 = fullfile(resdir,['compare_maxvalue_' responsespec1 '_' responsespec2 '_' tag '_boxstat.eps']);
saveas(H,filename);
set(H,'Renderer','painters')
saveas(H,filename2);
