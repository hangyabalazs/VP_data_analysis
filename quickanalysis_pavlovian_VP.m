function quickanalysis_pavlovian_VP(cellids, resdir_eg)
%QUICKANALYSIS2   Analysis of tetrode data.
%   QUICKANALYSIS2 is designed as an offline analysis tool for tetrode data
%   and behavior, which can be executed in an unsupervised manner on a
%   daily bases. It gives a quick overview of the experiment including
%   response profiles of clustered neurons, light-triggered PSTH and
%   psychometric plots of behavior. It relies on CellBase data handling
%   system.
%
%   QUICKANALYSIS2(ANIMALNO,SESSIONID,SESSIONSPEC) performs the analysis
%   for a session specified by the first two input arguments. SESSIONSPEC
%   should be a 1x3 logical array indicating the need for behavior,
%   recording and stimulation analysis.
%
%   QUICKANALYSIS2(ANIMALNO,SESSIONID,SESSIONSPEC,PROTOCOLTAG) accepts a
%   PROTOCOLTAG argument to allow calls to trial event conversion programs
%   for different versions of the behavioral protocol.
%
%   See also ADDNEWCELLS, PREALIGNSPIKES and VIEWCELL2B.

% Input argument check
narginchk(0,4)

% Animal, session
if isempty(cellids)
    cellids = { 'HDB13_170228a_2.1';
        'HDB13_170301a_3.2';
        'HDB32_181121a_5.3';
        'HDB32_181130a_2.2';
        'HDB32_181130a_4.1';
        'HDB25_180401a_7.3';
        };
end


% Stop if error
dbstop if error

% Directories
resdir = fullfile(resdir_eg,'example_neurons');   % results directory
if ~isdir(resdir)
    mkdir(resdir)
end



% Prealign spikes for trial events
problem_behav_cellid = [];
for iC = 1:length(cellids)
    cellid = cellids(iC);
    try
        prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_pavlovian,'filetype','event','ifsave',1,'ifappend',0)
    catch
        disp('Error in prealignSpikes.');
        problem_behav_cellid = [problem_behav_cellid cellid];
    end
end

%Excited by cue
G = figure;
pause(0.01)
viewcell2b(cellids(1),'TriggerName','StimulusOn','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'DeliverAllFeedback'}},'Partitions','#TrialType','window',[-3 3])
maximize_figure(G)

cellidt = char(cellids(1));
cellidt(cellidt=='.') = '_';
fnm = fullfile(resdir,[cellidt '_stimuluson_excitation.jpg']);   % save
saveas(G,fnm)
close(G)

%Inhibited by cue
G = figure;
pause(0.01)
viewcell2b(cellids(3),'TriggerName','StimulusOn','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'DeliverAllFeedback'}},'Partitions','#TrialType','window',[-3 3])
maximize_figure(G)

cellidt = char(cellids(3));
cellidt(cellidt=='.') = '_';
fnm = fullfile(resdir,[cellidt '_stimuluson_inhibition.jpg']);   % save
saveas(G,fnm)
close(G)

%Activated by reward
G = figure;
pause(0.01)
viewcell2b(cellids(5),'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#RewardedTrials','window',[-3 3])
maximize_figure(G)

cellidt = char(cellids(5));
cellidt(cellidt=='.') = '_';
fnm = fullfile(resdir,[cellidt '_reward_excitation.jpg']);   % save
saveas(G,fnm)
close(G)


%Inhibition by reward
G = figure;
pause(0.01)
viewcell2b(cellids(2),'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#RewardedTrials','window',[-3 3])
maximize_figure(G)

cellidt = char(cellids(2));
cellidt(cellidt=='.') = '_';
fnm = fullfile(resdir,[cellidt '_reward_inhibition.jpg']);   % save
saveas(G,fnm)
close(G)


G = figure;
pause(0.01)
viewcell2b(cellids(4),'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#PunishedTrials','window',[-3 3])
maximize_figure(G)

cellidt = char(cellids(4));
cellidt(cellidt=='.') = '_';
fnm = fullfile(resdir,[cellidt '_punishment_excitation.jpg']);   % save
saveas(G,fnm)
close(G)



G = figure;
pause(0.01)
viewcell2b(cellids(6),'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#PunishedTrials','window',[-3 3])
maximize_figure(G)

cellidt = char(cellids(6));
cellidt(cellidt=='.') = '_';
fnm = fullfile(resdir,[cellidt '_punishment_inhibition.jpg']);   % save
saveas(G,fnm)
close(G)