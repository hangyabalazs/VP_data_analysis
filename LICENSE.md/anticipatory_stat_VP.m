function anticipatory_stat_VP(animalID, resdir)
%ANTICIPATORY_STAT_VP compares anticipatory licking to two reinforcement
%predicting stimuli.
%   ANTICIPATORY_STAT_VP(ANIMALID, RESDIR) exctracts anticipatory licks in
%   a training session of an animal (ANIMALID) from TrialEvents (licks that
%   occur after the cue presentation and before reinforcement) and performs
%   a Wilcoxon signed rank test to compare the frequency of anticipatory
%   licking after reward or punishment predicting cue in Pavlovian
%   conditioning. The frist 20 trials can be free access of water in the
%   training protocol, therefore these trials are skipped. The results are
%   saved as an .xlsx file into RESDIR.
%
%   See also SIGNRANK.

%   Panna Hegedus, 05-Feb-2020
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu 

% Code review: BH 4/8/20

% Sessions
if nargin < 1
    mousepath = fullfile(getpref('cellbase','datapath'),'HDB18\');
else
    mousepath = fullfile(getpref('cellbase','datapath'),animalID, '\');
end
dr = dir(mousepath);
dr = dr(3:end);

% Make results directory
if ~isfolder(resdir)
    mkdir(resdir);
end

% Compare number of anticipatory licks for different trial types 
NumSessions = length(dr);
T = nan(NumSessions, 2);
for iS = 1:NumSessions  % loop through sessions
    sessionID = dr(iS).name;
    [TrialType1, TrialType2] = countlicks(animalID, sessionID);   % licks per trial type for each trial
    T(iS,1) = nansum(TrialType1)/(sum(~isnan(TrialType1))*1.2);   
    T(iS,2) = nansum(TrialType2)/(sum(~isnan(TrialType2))*1.2);
end
[p, h] = signrank(T(:,1), T(:,2));  %  perform Wilcoxon signed rank test across sessions
Results = [T; NaN NaN; p h ];  % write the result into an Excel file
xlswrite(fullfile(resdir, [animalID '_' sessionID '_behavior_stats.xlsx']),Results);
end

% -------------------------------------------------------------------------
function [g1, g2] = countlicks(animalID, sessionID)

% Load behavior data
cbdir = getpref('cellbase','datapath');
datapath = fullfile(cbdir,animalID,sessionID,'TrialEvents.mat');
VE = load(datapath);

if length(unique(VE.TrialType)) == 2
    ntrials = length(VE.NTrials); % length of each group
    
    % the frist 20 trials can be free access of water in the
    % training protocol, therefore these trials are skipped
    [g1, g2] = deal(zeros(1,ntrials));   % preallocate space
    for cT = 21:ntrials  % loop through trials from (skip first 20 trials)
        ws = VE.StimulusOn(cT);
        we = ws + 1.2; % window: stimulus duration + 200ms delay
        if VE.TrialType(cT) ==1 && (~isempty(VE.LickIn{cT})) % for TrialType = 1
            licks = VE.LickIn{cT}; % lick
            licks = licks(licks>ws);
            licks = licks(licks<we);
            g1(1, cT) = length(licks);
        elseif VE.TrialType(cT) ==2 && (~isempty(VE.LickIn(cT))) % for TrialType = 2
            licks = VE.LickIn{cT};
            licks = licks(licks>ws);
            licks = licks(licks<we);
            g2(1, cT) = length(licks);
        end
    end
else
    g1 = NaN;
    g2 = NaN; % do not make statistics, if there is only 1 trialtype
end
end