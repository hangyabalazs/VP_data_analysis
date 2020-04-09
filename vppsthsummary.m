function vppsthsummary
%VPPSTHSUMMARY   PSTH summary figures.
%   VBPSTHSUMMARY plots population PSTH for VP neurons recorded in a
%   pavlovian task. Z-scored PSTHs are used throughout the analysis.
%   Average PSTHs aligned to behavioral events are plotted.
%
%   See also NBPSTHSUMMARY.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   8-March-2014

% Directories
global DATAPATH
resdir = ([DATAPATH 'psthsummary\punish\']);   % result directory
if ~isdir(resdir)
    mkdir(resdir)
end

% Cells
selstr = '"ispv"==1';
cellids = selectcell(selstr);   % select based on responses to behav. events
NumCells = length(cellids);

% LED-on PSTH
allpsth = [];
allstats = [];
wn = [-0.5 1];
dt = 0.001;
time = wn(1):dt:wn(2);   % time vector
for k = 1:NumCells
    cellid = cellids{k};   % cell ID
    disp(cellid)
    
    % Align event
    alignevent = 'DeliverAllFeedback';
    alignfilter = 'Punishment==2';
    
    % Calcualte PSTH
    [psth, spsth, ~, ~, spt] = ultimate_psth(cellid,'trial',...
        alignevent,wn,'event_filter','custom','filterinput',alignfilter,...
        'dt',dt,'display',true,'sigma',0.02,'parts','all','isadaptive',2,...
        'maxtrialno',Inf);
%     psth = viewcell2b(cellid,'TriggerName','LEDOn','SortEvent','TrialStart','eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','all','window',wn);
    
    % Concatenate data from different cells
    allpsth = [allpsth; spsth]; %#ok<AGROW>
end

% Normalize
allpsth = zscore(allpsth,0,2);
allpsth = zscore(allpsth,0,2);

% Population PSTHs for groups
figure   % images proportional to number of cells in the groups
imagesc(allpsth);
colormap hot
set(gca,'CLim',[-1 20])
saveas(gcf,fullfile(resdir,'poppsth_PunishmentOn.fig'))

[m1 m2] = max(allpsth,[],2);
[srt Ia] = sort(m2,'ascend');   % sort based on response latency
figure   % plot all PSTHs, sorted
imagesc(time,1:size(allpsth,1),allpsth(Ia,:))
colormap(hot)
saveas(gcf,fullfile(resdir,'poppsth_PunishmentOn_sorted.fig'))

% Average PSTH
figure
hold on;
baseline = mean(allpsth(:,1:20));
mn = mean(baseline);   % baseline for feedback alignment
errorshade(time,mean(allpsth)-mn,std(allpsth)/sqrt(size(allpsth,1)),...
    'LineColor','k','ShadeColor','k')

ymx = 6;
set(gca,'YLim',[-2 ymx],'YTick',[0 ymx/2 ymx],'YTickLabel',{'0' '' num2str(ymx)});
line([0 0],ylim,'Color','k','LineStyle',':')
line(xlim,[0 0],'Color','k')
xlabel('Time from punishment (s)')
ylabel({'Normalized';'firing rate'})
setmyplot_balazs
saveas(gcf,fullfile(resdir,'psth_average_PunishmentOn2.fig'))
set(gcf,'Renderer','painters')
saveas(gcf,fullfile(resdir,'psth_average_PunishmentOn2.eps'))