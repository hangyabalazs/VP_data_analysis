function avg_psth_VP(cellids1,cellids2,responsespec,resdir,spike)
%AVG_PSTH_VP   Average normalized peri-event time histograms.
%   AVG_PSTH_VP(CELLIDS1,CELLIDS2,RESPONSESPEC) calculates normalized
%   average PSTHs aligned to RESPONSESPEC (cue, reward or punishment) for
%   CELLIDS1 and CELLIDS2, respectively. Result plots are generated by
%   ULTIMATE_PSTH (CellBase) and saved in RESDIR. Indivual PSTHs are
%   z-scored by the mean and SD of a predefined baseline window before
%   averaging.
%
%   AVG_PSTH_VP(CELLIDS,CELLIDS,RESPONSESPEC,'BURST') compares the PSTHs of
%   bursts and single spikes.
%
%   See also ULTIMATE_PSTH.

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   03-Dec-2019

%   Code review: BH 12/6/19, 4/8/20, 10/16/20

% Directories
if ~isfolder(resdir)
    mkdir(resdir)
end

% Input argument check
if nargin < 5
    spike_def = {[] []};
elseif strcmp(spike,'burst')
    spike_def = {'burst' 'single'};
end

% Aligning: cue, reward or punishment
switch responsespec
    case {'cue','cueresponse'}
        alignevent = 'StimulusOn';
        filterevent = 'StimulusOn';
        bwin = [-1 0]; % baseline window
    case {'reward','rewardresponse'}
        alignevent = 'DeliverAllFeedback';
        filterevent = 'RewardedTrials==1';
        bwin = [-3 -2]; % baseline window
    case {'punishment','punishresponse','punishmentresponse'}
        alignevent = 'DeliverAllFeedback';
        filterevent = 'PunishedTrials==1';
        bwin = [-3 -2]; % baseline window
    case {'omission','omissionresponse'}
        alignevent = 'DeliverAllFeedback';
        filterevent = 'Omission==1|Omission==2';
        bwin = [-3 -2]; % baseline window
    otherwise
        error('MATLAB:avg_psth_VP:responsespecUnknown','Unrecognized input value for responsespec.')
end
disp(['Aligning to: ' alignevent])

% Time vector
wn = [-4 4]; % time window
dt = 0.001; % time resolution
time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector
btime = bwin(1)*1000:dt*1000:bwin(2)*1000;   % baseline time vector
[~, baseline_inx] = intersect(time,btime);   % baseline indices

% Ultimate PSTH
R1 = runanalysis(@ultimate_psth,...
    'trial', alignevent, wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',0,...
    'event_filter','custom', 'filterinput',filterevent,'spike_def', spike_def{1}, 'maxtrialno',Inf,...
    'baselinewin',bwin,'testwin',[0 0.5],'relative_threshold',0.1,'cellids',cellids1);
R2 = runanalysis(@ultimate_psth,...
    'trial', alignevent, wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',0,...
    'event_filter','custom', 'filterinput',filterevent,'spike_def', spike_def{2}, 'maxtrialno',Inf,...
    'baselinewin',bwin,'testwin',[0 0.5],'relative_threshold',0.1,'cellids',cellids2);
R = [{R1} {R2}];
Rinx = 1;

% Normalization
r1 = cell2mat(R1(:,Rinx));
bl = r1(:,baseline_inx);   % baseline matrix
psth_R1 = (r1 - repmat(mean(bl,2),1,size(r1,2))) ./ repmat(std(bl,[],2),1,size(r1,2));
r2 = cell2mat(R2(:,Rinx));
bl = r2(:,baseline_inx);   % baseline matrix
psth_R2 = (r2 - repmat(mean(bl,2),1,size(r2,2))) ./ repmat(std(bl,[],2),1,size(r2,2));

% Plot color coded PSTH
Ha = figure;
[~, sinx] = sort(max(psth_R1,[],2), 'descend');   % sort by maximum
imagesc(psth_R1(sinx,:))
colorbar
Hb = figure;
[~, sinx] = sort(min(psth_R2,[],2), 'descend');   % sort by minimum
imagesc(psth_R2(sinx,:))
colorbar

H = figure;
errorshade(time,nanmean(psth_R1),nanstd(psth_R1)/sqrt(size(psth_R1,1)),...
    'LineColor',[0 0 1],'ShadeColor',[0 0 1]) % excitation
errorshade(time,mean(psth_R2),std(psth_R2)/sqrt(size(psth_R2,1)),...
    'LineColor',[0 0 0.5],'ShadeColor',[0 0 0.5]) % inhibition
set(gca,'XLim',[-3000 3000]);
set(gca,'YLim',[-10 20]);
set(gcf, 'renderer', 'painters')

% Save
fnm = fullfile(resdir,[responsespec '_average.jpg']);   % jpg
fnm2 = fullfile(resdir,[responsespec '_average.fig']);   % fig
fnm3 = fullfile(resdir,[responsespec '_average.mat']);   % full PSTH matrix
saveas(H,fnm);
saveas(H,fnm2);
save(fnm3,'R');
close(H)

fnma = fullfile(resdir,[responsespec '_image1.jpg']);   % sorted normalized PSTHs
fnmb = fullfile(resdir,[responsespec '_image2.jpg']);   % sorted normalized PSTHs
saveas(Ha,fnma);
saveas(Hb,fnmb);
close(Ha,Hb)