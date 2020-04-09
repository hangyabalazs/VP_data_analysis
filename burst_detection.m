function burst_detection(cellids, spiketype)
%BURST_DETECTION finds bursts within a spike train.
%
%   BURST_DETECTION(CELLIDS,SPIKETYPE )detects bursts (SPIKETYPE = 'burst') or
%   single spikes (SPIKETYPE = 'single') of cells defined by CELLID. The
%   PSTHs of bursts or single spikes are generated with VIEWCELL2BURSTS
%   or VIEWCELL2P respectively. Event aligned bursts or single spikes are also
%   saved into the corresponding folder.
%
%   See also REALIGNBURSTS, PREALIGNSINGLESPIKES, VIEWCELL2BURST and
%   VIEWCELL2P.

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   05-Feb-2020

%   Code review: BH 2/11/20

% Input argument check
narginchk(0,2);
if nargin < 1  %if no cells pre-selected
    vpcells = vpselectcells([]);
else
    vpcells = cellids;
end

% Directories
resdir = fullfile(getpref('cellbase','datapath'),'VP','burst_psth2');   % results directory
if ~isdir(resdir)
    mkdir(resdir)
end

% Pre-align bursts and single spikes
numCells = length(vpcells);   % number of VP neurons
switch spiketype
    case 'burst'
        for i = 1:numCells
            
            problem_stim_cellid = [];
            cellid = vpcells(i);
            try   % prealign bursts
                prealignBursts(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_pavlovian,'filetype','event','ifsave',1,'ifappend',0)
            catch
                disp('Error in prealignBursts.');
                problem_stim_cellid = [problem_stim_cellid cellid];
            end
        end
        
        % Plot and save
        for k = 1:numCells
            G = figure;
            pause(0.01)
            viewcell2burst(vpcells(k),'TriggerName','StimulusOn','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'DeliverAllFeedback'}},'Partitions','#TrialType','window',[-5 5])
            maximize_figure(G)
            
            cellidt = vpcells{k};
            cellidt(cellidt=='.') = '_';
            fnm = fullfile(resdir,[cellidt '_stimuluson_burst.jpg']);   % save
            saveas(G,fnm)
            set(G,'Renderer','painters');
            saveas(G,fullfile(resdir,[cellidt '_stimuluson_burst.eps']));
            close(G)
        end
        
        for k = 1:numCells
            J = figure;
            pause(0.01)
            viewcell2burst(vpcells(k),'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#TrialType','window',[-5 5]);
            maximize_figure(J)
            cellidt = vpcells{k};
            cellidt(cellidt=='.') = '_';
            fnm = fullfile(resdir,[cellidt '_DeliverAllFeedback_burst.jpg']);   % save
            saveas(J,fnm)
            set(J,'Renderer','painters');
            saveas(J,fullfile(resdir,[cellidt '_DeliverAllFeedback_burst.eps']));
            close(J)
        end
        
        % Test activation/inhibition for cue, reward and punishment
        vpresponsesorter_burst(vpcells,1,resdir, 'cue');
        vpresponsesorter_burst(vpcells,1,resdir, 'rew');
        vpresponsesorter_burst(vpcells,1,resdir, 'pun');
        
    otherwise 'single'
        for i = 1:numCells
            
            problem_stim_cellid = [];
            cellid = vpcells(i);
            try   % prealign single spikes
                prealignSinglespikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_pavlovian,'filetype','event','ifsave',1,'ifappend',0)
            catch
                disp('Error in prealignBursts.');
                problem_stim_cellid = [problem_stim_cellid cellid];
            end
        end
        
        % Plot and save
        for k = 1:numCells
            G = figure;
            pause(0.01)
            viewcell2p(vpcells(k),'TriggerName','StimulusOn','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'DeliverAllFeedback'}},'Partitions','#TrialType','window',[-5 5])
            maximize_figure(G)
            
            cellidt = vpcells{k};
            cellidt(cellidt=='.') = '_';
            fnm = fullfile(resdir,[cellidt '_stimuluson_single.jpg']);   % save
            saveas(G,fnm)
            set(G,'Renderer','painters');
            saveas(G,fullfile(resdir,[cellidt '_stimuluson_single.eps']));
            close(G)
        end
        
        for k = 1:numCells
            J = figure;
            pause(0.01)
            viewcell2p(vpcells(k),'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#TrialType','window',[-5 5]);
            maximize_figure(J)
            cellidt = vpcells{k};
            cellidt(cellidt=='.') = '_';
            fnm = fullfile(resdir,[cellidt '_DeliverAllFeedback_single.jpg']);   % save
            saveas(J,fnm)
            set(J,'Renderer','painters');
            saveas(J,fullfile(resdir,[cellidt '_DeliverAllFeedback_single.eps']));
            close(J)
        end
        
        % Test activation/inhibition for cue, reward and punishment
        vpresponsesorter_single(vpcells,1,resdir, 'cue'); 
        vpresponsesorter_single(vpcells,1,resdir, 'rew');
        vpresponsesorter_single(vpcells,1,resdir, 'pun');
end