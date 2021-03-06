function avg_psth_VP_animals(cellids, isbehav, responsespec,resdir)
%AVG_PSTH_VP_ANIMALS   Average normalized peri-event time histograms.
%   [MAXVAL, MINVAL] = AVG_PSTH_VP_ANIMALS(CELLIDS,RESPONSESPEC, RESDIR)
%   calculates normalized average PSTHs aligned to RESPONSESPEC (cue,
%   reward or punishment) for each animal, respectively. If ISBEHAV is
%   true, only animals with good behavioral measurements are analyzed.
%   MAXVAL and MINVAL (maximal and minimal neuronal response) are
%   calculated by ULTIMATE_PSTH. Result plots are also generated by
%   ULTIMATE_PSTH (CellBase) and saved in RESDIR. Indivual PSTHs are
%   z-scored by the mean and SD of a predefined baseline window before
%   averaging.
%
%   See also ULTIMATE_PSTH.

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   16-June-2020

% Directories
if ~isfolder(resdir)
    mkdir(resdir)
end

% Aligning: cue, reward or punishment
switch responsespec
    case {'cue'}
        alignevent = 'StimulusOn';
        filterevent = 'StimulusOn';
        bwin = [-1 0]; % baseline window
        twin = [0 0.5];
    case {'reward'}
        alignevent = 'DeliverAllFeedback';
        filterevent = 'RewardedTrials==1';
        bwin = [-3 -2]; % baseline window
        twin = [0 0.2];
    case {'punishment'}
        alignevent = 'DeliverAllFeedback';
        filterevent = 'PunishedTrials==1';
        bwin = [-3 -2]; % baseline window
        twin = [0 0.2];
    otherwise
        error('MATLAB:avg_psth_VP:responsespecUnknown','Unrecognized input value for responsespec.')
end
disp(['Aligning to: ' alignevent])

% Select neurons
if isbehav
    animals = {'HDB13' 'HDB32'}
    cellids = rmfield(cellids, {'HDB25' 'HDB17'});
    cell_groups = struct2cell(cellids);
    labels1 = 'good behav animals'; % legend for later plotting
    labels2 = 'good behav animals'; % legend for later plotting
else
    animals = fieldnames(cellids); % animal name
    cell_groups = struct2cell(cellids);
    labels1 = animals; % legend for later plotting
    labels2 = animals; % legend for later plotting
end

% Time vector
wn = [-4 4]; % time window
dt = 0.001; % time resolution
time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector
btime = bwin(1)*1000:dt*1000:bwin(2)*1000;   % baseline time vector
[~, baseline_inx] = intersect(time,btime);   % baseline indices

% Ultimate PSTH
for i = 1:length(animals)
    aID = animals{i};
    cellids1 = cell_groups{i}.(responsespec).excitation;
    cellids2 = cell_groups{i}.(responsespec).inhibition;
    R1.(aID) = runanalysis(@ultimate_psth,...
        'trial', alignevent, wn,...
        'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
        'event_filter','custom', 'filterinput',filterevent,'maxtrialno',Inf,...
        'baselinewin',bwin,'testwin',twin, 'relative_threshold',0.1,'cellids',cellids1);
    R2.(aID) = runanalysis(@ultimate_psth,...
        'trial', alignevent, wn,...
        'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
        'event_filter','custom', 'filterinput',filterevent,'maxtrialno',Inf,...
        'baselinewin',bwin,'testwin',twin,'relative_threshold',0.1,'cellids',cellids2);
end

Rinx = 2; % choose spsth

% Normalization
if isbehav
    [psth_R1, psth_R2, exc_stat, inh_stat] = deal(cell(1,1)); % prealloace space
    
    if isempty(R1.HDB13) % if one of the animals does not have cells with specified response
        all_cells_behav1 = [R1.HDB32(:,Rinx)];
    elseif isempty (R1.HDB32)
        all_cells_behav1 = [R1.HDB13(:,Rinx)];
    elseif ~isempty(R1.HDB13) && ~isempty(R1.HDB32)
        all_cells_behav1 = [R1.HDB13(:,Rinx) ; R1.HDB32(:,Rinx)];
    end
    
    r1 = cell2mat(all_cells_behav1);
    bl = r1(:,baseline_inx);   % baseline matrix
    psth_R1{1} = (r1 - repmat(mean(bl,2),1,size(r1,2))) ./ repmat(std(bl,[],2),1,size(r1,2));
    
    if isempty (R2.HDB13) % if one of the animals does not have cells with specified response
        all_cells_behav2 = [R2.HDB32(:,Rinx)];
    elseif isempty (R2.HDB32)
        all_cells_behav2 = [R2.HDB13(:,Rinx)];
    elseif  ~isempty(R2.HDB13) && ~isempty(R2.HDB32)
        all_cells_behav2 = [R2.HDB13(:,Rinx) ; R2.HDB32(:,Rinx)];
    end
    
    r2 = cell2mat(all_cells_behav2);
    bl = r2(:,baseline_inx);   % baseline matrix
    psth_R2{1} = (r2 - repmat(mean(bl,2),1,size(r2,2))) ./ repmat(std(bl,[],2),1,size(r2,2));
else
    [psth_R1, psth_R2, exc_stat, inh_stat] = deal(cell(1,length(animals)));
    for k = 1:length(animals) %excitation - average PSTH
        if ~isempty(R1.(animals{k}))
            r1 = cell2mat(R1.(animals{k})(:,Rinx));
            bl = r1(:,baseline_inx);   % baseline matrix
            psth_R1{k} = (r1 - repmat(mean(bl,2),1,size(r1,2))) ./ repmat(std(bl,[],2),1,size(r1,2));
            exc_stat{k} = R1.(animals{k})(:,6);
        end
    end
    
    for k = 1:length(animals) % inhibition - average PSTH
        if ~isempty(R2.(animals{k}))
            r2 = cell2mat(R2.(animals{k})(:,Rinx));
            bl = r2(:,baseline_inx);   % baseline matrix
            psth_R2{k} = (r2 - repmat(mean(bl,2),1,size(r2,2))) ./ repmat(std(bl,[],2),1,size(r2,2));
            inh_stat = R2.(animals{k})(:,6);
        end
    end
%     
%     exc_stat = exc_stat(~cellfun(@isempty, psth_R1));
%     inh_stat = inh_stat(~cellfun(@isempty, psth_R2));

    labels1 = labels1(~cellfun(@isempty, psth_R1));
    labels2 = labels2(~cellfun(@isempty, psth_R2));
    
    psth_R1 = psth_R1(~cellfun(@isempty, psth_R1)); % remove empty cells
    psth_R2 = psth_R2(~cellfun(@isempty, psth_R2));
end


H = figure;
colors = {[1 0 0] [0 1 0] [0 0 1] [1 0 1]};
for j = 1:length(psth_R1)
    errorshade(time,mean(psth_R1{j}),std(psth_R1{j})/sqrt(size(psth_R1{j},1)),...
        'LineColor',colors{j},'ShadeColor',colors{j}) % excitation
    hold on
end
legend(labels1);
set(gca,'XLim',[-3000 3000]);
set(gca,'YLim',[-5 50]);

G = figure;
for m = 1:length(psth_R2)
    errorshade(time,mean(psth_R2{m}),std(psth_R2{m})/sqrt(size(psth_R2{m},1)),...
        'LineColor',colors{m},'ShadeColor',colors{m}) % inhibition
end
legend(labels2);
set(gca,'XLim',[-3000 3000]);
set(gca,'YLim',[-20 5]);

% Save
fnm = fullfile(resdir,[responsespec 'excitation_average.jpg']);   % jpg
fnm2 = fullfile(resdir,[responsespec 'excitation_average.fig']);   % fig
fnm3 = fullfile(resdir,[responsespec 'excitation_average.mat']);   % full PSTH matrix
saveas(H,fnm);
saveas(H,fnm2);
save(fnm3,'R1');
close(H)

fnm = fullfile(resdir,[responsespec 'inhibition_average.jpg']);   % jpg
fnm2 = fullfile(resdir,[responsespec 'inhibition_average.fig']);   % fig
fnm3 = fullfile(resdir,[responsespec 'inhibition_average.mat']);   % full PSTH matrix
saveas(G,fnm);
saveas(G,fnm2);
save(fnm3,'R2');
close(G)