function lat = latency_dist(vpcells,sourcedir)
%LATENCY_DIST   Load response latency statistics.
%   LAT = LATENCY_DIST(CELLIDS,DIR) loads latency statistics for cue,
%   reward and punishment responsive cells of CELLIDS from SOURCEDIR. The results
%   are returned in a struct containing one field for each of cue,
%   reward and punishment. These individual fields contain an array of
%   five fields for start, peak, duration, maxvalue and minvalue for activated and inhibited
%   neurons, respectively.
%
%   See also ULTIMATE_PSTH and PSTH_STATS.

%   Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hegedus.panna@koki.mta.hu
%   19-Apr-2019

%   Code review: BH 12/20/19, 1/3/20

% Input arguments
if nargin < 1 || isempty(vpcells)
    vpcells = vpselectcells(vpcells);
end

% Responsive cells
cueresp = getvalue('cueresponse', vpcells);
rewardresp = getvalue('rewardresponse', vpcells);
punishmentresp = getvalue('punishresponse', vpcells);

% Load latency statistics
stimuli = [{'cue'}, {'reward'}, {'punishment'}];
lat = struct;
for i = 1:length(stimuli)  % loop through cue, reward and punishment
    s = stimuli{i};
    
    switch s
        case 'cue'
            p1_vpcells = vpcells(cueresp == 1); % excited by cue
            p2_vpcells = vpcells(cueresp == -1); % inhibited by cue
            stimstr = '_StimulusOn_TrialType.mat';
            
        case 'reward'
            p1_vpcells = vpcells(rewardresp == 1); % excited by reward
            p2_vpcells = vpcells(rewardresp == -1); % inhibited by reward
            stimstr = '_DeliverAllFeedback_RewardedTrials.mat';
            
        case 'punishment'
            p1_vpcells = vpcells(punishmentresp == 1); % excited by punishment
            p2_vpcells = vpcells(punishmentresp == -1); % inhibited by punishment
            stimstr = '_DeliverAllFeedback_PunishedTrials.mat';
    end
    
    response = [{'exc'},{'inh'}];  % response types
    cells = [{p1_vpcells}, {p2_vpcells}];
    numResponses = length(response);
    [latency_start, latency_peak, latency_duration, minvalue, maxvalue] = deal(cell(1,numResponses));
    for j = 1:numResponses  % loop through excitation and inhibition
        r = response{j};
        
        currentCells = cells{j};
        numCells = length(currentCells);
        [latency_start{j}, latency_peak{j}, latency_duration{j}] = deal(nan(1,numCells));
        for k = 1:numCells   % loop through cells
            cellid = char(currentCells(k));
            [animalID, sessionID, tt, u]  = cellid2tags(cellid);
            sttc = [animalID '_' sessionID '_' num2str(tt) '_' num2str(u)];
            load([sourcedir '\' sttc stimstr]);   % load response analysis file
            close all;
            
            if strcmp(r,'exc') == 1 % if excited
                latency_start{j}(k) = stats1{1,1}.activation_start(1); %#ok<USENS>
                latency_peak{j}(k) = stats1{1,1}.activation_peak(1);
                latency_duration{j}(k) = stats1{1,1}.activation_time(1);
                maxvalue{j}(k) = stats1{1,1}.maxvalue(1);
                minvalue{j}(k) = stats1{1,1}.minvalue(1);
            elseif strcmp(r,'inh') == 1 % if inhibited
                latency_start{j}(k) = stats1{1,1}.inhibition_start(1);
                latency_peak{j}(k) = stats1{1,1}.inhibition_peak(1);
                latency_duration{j}(k) = stats1{1,1}.inhibition_time(1);
                maxvalue{j}(k) = stats1{1,1}.maxvalue(1);
                minvalue{j}(k) = stats1{1,1}.minvalue(1);
            end
        end
    end
    
    % Output
    lat.(s).excited.latency_start = latency_start{1};
    lat.(s).inhibited.latency_start = latency_start{2};
    lat.(s).excited.latency_peak = latency_peak{1};
    lat.(s).inhibited.latency_peak = latency_peak{2};
    lat.(s).excited.latency_duration = latency_duration{1};
    lat.(s).inhibited.latency_duration = latency_duration{2};
    lat.(s).excited.maxvalue = maxvalue{1};
    lat.(s).inhibited.maxvalue = maxvalue{2};
    lat.(s).excited.minvalue = minvalue{1};
    lat.(s).inhibited.minvalue = minvalue{2}; 
end