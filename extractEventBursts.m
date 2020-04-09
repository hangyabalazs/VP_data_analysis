function [event_stimes,event_windows] = extractEventBursts(cellid,events,type,TE)
%EXTRACTEVENTBURSTS   Extract and align singe spikes or bursts in each trial relative to trial events.
%   [EVENT_STIMES,EVENT_WINDOWS] = EXTRACTEVENTBURSTS(CELLID,EVENTS,TYPE)
%   Input arguments:
%       CELLID  - ID of a unit in CellBase
%       EVENTS  - Nx4 cell array of events (see DEFINEEVENTSEPOCHS)
%       TPYE - 'burst' extracts bursts from the spike trian, 'single'
%           extracts single spikes
%
%   Output arguments:
%       EVENT_STIMES - 1xN cell array of 1xTrialNum cell arrays of spike times
%           relative to the first event trigger in 'events'.
%       EVENT_WINDOW - 1xN cell array of 1xTrialNum arrays noting the time window 
%           of extraction for each trial.
%
%   See also EXTRACTEVENTBURSTS and PREALIGNSPIKES.

%   Edit log: AK 11/06, ZFM 10/04, BH 6/27/11, BH 2/11/20

% Initialize
disp(['Extracting event-triggered spike times for ' char(cellid)])
NUMevents = size(events,1);
event_stimes = cell(NUMevents,1);
event_windows = cell(NUMevents,1);

% Load TrialEvents
if nargin < 4
    TE = loadcb(cellid,'Events');   
end

% Load spike times
stimes = loadcb(cellid,'Spikes');
% stimes = stimes * 1e-4;   % conversion factor into seconds

% Inter-spike intervals
isi = diff(stimes);  % inter-spike interval
    burstinx = find(isi<0.01);   % ISI < 10 ms
    
    % Detect bursts
    % Burst: first ISI < 10 ms, subsequent ISIs < 15 ms
    bursts = {};
    used = [];
    burst1 = {};
    for k = 1:length(burstinx)
        if ismember(burstinx(k),used)
            continue   % already included in the previous burst
        end
        bursts{end+1} = stimes([burstinx(k), burstinx(k)+1]); 
        next = burstinx(k) + 1;
        if next > length(isi)
            burst1{end+1} = bursts{end}(1);
            continue   % last ISI of the cell
        end
        while isi(next) < 0.015  % if conseq. ISIs < 15 ms
            bursts{end} = [bursts{end}; stimes(next+1)];
            used = [used next]; %#ok<AGROW>
            next = next + 1;
            if next > length(isi)
                break   % last ISI of the cell, quit while loop
            end
        end
        burst1{end+1} = bursts{end}(1);
    end
    bursts = vertcat(bursts{:});
    single = setdiff(stimes, bursts);
    burst1 = cell2mat(burst1)';


switch type
    case 'single'
        spike_times = single;
    case 'burst'
        spike_times = bursts;
end

for iE = 1:NUMevents
    
    % Event references 
    TriggerEvent1 = events{iE,2}; 
    TriggerEvent2 = events{iE,3};
    
    % Window 
    offset1 = events{iE,4}(1);
    offset2 = events{iE,4}(2);
    
    if isfield(TE,TriggerEvent1) && isfield(TE,TriggerEvent2)
        NUMtrials = length(TE.TrialStart);
        if ~iscell(TE.(TriggerEvent1))    % original code
            EventTimes1 = TE.(TriggerEvent1) + TE.TrialStart + offset1;
            EventTimes2 = TE.(TriggerEvent2) + TE.TrialStart + offset2;
            
            % Note event windows are relative to TriggerEvent1
            event_windows{iE}(1,1:NUMtrials) = offset1;
            event_windows{iE}(2,:) = (EventTimes2 - EventTimes1) + offset1;
            
            trial_stimes = cell(1,NUMtrials);
            for iT = 1:NUMtrials
                if ~isnan(EventTimes1(iT)) && ~isnan(EventTimes2(iT))
                    trial_stimes{iT} = spike_times(spike_times>=EventTimes1(iT)&spike_times<EventTimes2(iT)) ...
                        - EventTimes1(iT) + offset1;
                end
            end   % iT
            event_stimes{iE} = trial_stimes;
            
        else  % addition by BH to handle lick events that are cells arrays with a cell of time stamps corresponding to each trial
            EventTimes1 = TE.(TriggerEvent1);
            EventTimes2 = TE.(TriggerEvent2);
            event_windows{iE} = cell(1,NUMtrials);
            trial_stimes = cell(1,NUMtrials);
            for iT = 1:NUMtrials
                NUMlicks = length(EventTimes1{iT});
                EventTimes1{iT} = EventTimes1{iT} + TE.TrialStart(iT) + offset1;
                EventTimes2{iT} = EventTimes2{iT} + TE.TrialStart(iT) + offset2;
                event_windows{iE}{iT}(1,1:NUMlicks) = offset1; % note event windows are relative to TriggerEvent1
                event_windows{iE}{iT}(2,:) = (EventTimes2{iT} - EventTimes1{iT}) + offset1;
                
                trial_stimes{iT} = cell(1,NUMlicks);
                for iL = 1:NUMlicks
                    trial_stimes{iT}{iL} = spike_times(spike_times>=EventTimes1{iT}(iL) & ...
                        spike_times<EventTimes2{iT}(iL))-EventTimes1{iT}(iL)+offset1;
                end
            end
            event_stimes{iE} = trial_stimes;
        end
        
    else % some events not found
        disp(iE)
        disp(sprintf('Warning: %s missing events: %s and/or %s',char(cellid),TriggerEvent1,TriggerEvent2));
        event_stimes{iE} = [];
        event_windows{iE} = [NaN NaN];
    end
end   % iE
