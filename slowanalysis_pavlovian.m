function slowanalysis_pavlovian(vpcells, preprocess, fig1spec, fig2spec, ...
    fig3spec, fig4spec, fig5spec, fig6spec, fig7spec, fig8spec, suppl_figspec)
%SLOWANALYSIS_PAVLOVIAN   Analysis main function for VP project.
%   SLOWANALYSIS_PAVLOVIAN(VPCELLS, PRE, F1, F2, ..., F8, S1) is the main
%   function for ventral pallidal data analysis. Data analysis is performed
%   for CellIDs in VPCELLS - or all well-isolated units determined by
%   VPSELECTCELLS if the first input is empty. The analysis modules are
%   controlled by logical input variables grouped by figures (and
%   pre-processing including response analysis, auto-correlation and
%   cross-correlation). This code is for further analysis after running
%   QUICKANALYSIS2_P on the CellBase dataset.
%
%   See also QUICKANALYSIS2_P, VPSELECTCELLS, VPRESPONSESORTER, VPACG and
%   VPCCG.

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   05-Feb-2020

%   Code review: BH 2/12/20, 4/8/20

% Input argument check
narginchk(0,11)
if nargin < 1
    vpcells = [];
end

% Choose CellBase
usr = getenv('username');
if ismember(usr,{'hangya.balazs','hangyab','Hangya Balázs'})
    choosecb('auditory_pavlovian_cellbase')   % choose CellBase
else
    choosecb('VP_CellBase')   % choose CellBase
end

% Control which analysis modules to perform
if nargin < 2
    preprocess = [1 1 1];
end
response_analysis = preprocess(1);   % RESPONSE SORTER
perform_acg = preprocess(2);   % ACG
perform_ccg = preprocess(3);   % CCG
if nargin < 3
    fig1spec = [1 1 1];
end
lick_example = fig1spec(1); % FIG1.E-F - example lick raster
lick_average = fig1spec(2); % FIG1.G - average lick raster
lick_stat = fig1spec(3); % FIG1.H - performance statistics per animal
if nargin < 4
    fig2spec = [1 1 1 1];
end
response_example = fig2spec(1); % FIG2.A-F - example PSTHs
response_pie = fig2spec(2); % FIG2.G-I - pie charts
response_PSTH = fig2spec(3); % FIG2.J-L, O - average PSTHs and scplot
box_stat = fig2spec(4); % FIG2.M-N - boxplots for latency peak
if nargin < 5
    fig3spec = [1 1];
end
compare_exp_psth = fig3spec(1); % FIG.3A,C,E - compare responses to expected and unexpected outcomes, PSTH
compare_exp_box = fig3spec(2); % FIG.3B,D,F - compare responses to expected and unexpected outcomes, boxstat
if nargin < 6
    fig4spec = [1 1];
end
autocorr_bursting = fig4spec(1); % FIG.4A - autocorrelation examples
pie_PSTH_bursting = fig4spec(2); % FIG.4B-I - pie charts and PSTHs
if nargin < 7
    fig5spec = [1 1];
end
autocorr_rhythmic = fig5spec(1); % FIG.5A - autocorrelation examples
pie_PSTH_rhythmic = fig5spec(2); % FIG.5B-J - pie charts and PSTHs
if nargin < 8
    fig6spec = [1 1 1];
end
single_spikes = fig6spec(1); % FIG.6A,C,D,F - raster plot and PSTH for single spikes
bursts = fig6spec(2); % FIG.6B,E - raster plot and PSTH for bursts
pie_burst_single = fig6spec(3); % FIG.6G-H
if nargin < 9
    fig7spec = 1;
end
select_pairs = fig7spec(1); % FIG.7A-B - network of synaptic pairs
if nargin < 10
    fig8spec = [1 1 1];
end

if nargin < 11  % do supplementary analysis or not
    suppl_figspec = [1 1 1 1 1];
end
figS1spec = suppl_figspec(1); % FigS1
figS2spec = suppl_figspec(2); % FigS2
figS3spec = suppl_figspec(3); % FigS3
figS4spec = suppl_figspec(4); % FigS4
figS5spec = suppl_figspec(5); % FigS5

vpcross = fig8spec(1); % FIG.8A-B - example and average CCG
sync_pie_PSTH = fig8spec(2); % FIG.8C-I - pie chart and average PSTH for sync and non-sync
norhythmandburst = fig8spec(3); % FIG.8K - pie chart of sync firing in the non-bursting-non-rhythmic group

% Stop if error
dbstop if error

% Directories
resdir = [getpref('cellbase','datapath') '\_paper_figs\code_review3\']; % directory for generated images
if ~isfolder(resdir)
    mkdir(resdir)
end

% Getting cellids for VP neurons
loadcb
vpcells = vpselectcells(vpcells);

% Grouping neurons based on their response during behavior - properties added to TheMatrix
response_resdir = fullfile(resdir,'vpresponsesorter');   % results directory
if response_analysis
    vpresponsesorter(vpcells,1,response_resdir,'cue');
    vpresponsesorter(vpcells,1,response_resdir,'rew');
    vpresponsesorter(vpcells,1,response_resdir,'pun');
end

% Auto-correlation
acg_resdir = fullfile(resdir,'vpacg');   % results directory
if perform_acg
    vpacg(vpcells,acg_resdir,true);
end

% Cross-correlation
ccg_resdir = fullfile(resdir,'vpccg');   % results directory
if perform_ccg
    vpccg(vpcells,ccg_resdir,true);
end

% Pre-processing directories
responsedir = fullfile(getpref('cellbase','datapath'),'VP','vpresponsesorter5');  % 'response_resdir'
acgdir = acg_resdir;
ccgdir = ccg_resdir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example lick rasters - FIG.1E-F
if lick_example
    figure;  % example lick PSTH + raster
    viewlick({'HDB13' '170316a'},'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav',...
        'ShowEvents',{{'StimulusOn' 'StimulusOff'}},...
        'Partitions','#TrialType','window',[-5 5]   )% lick rasters aligned to stimulus onset
end

% Average lick rasters - FIG.1G
if lick_average
    resdir_lick_psth = fullfile(resdir,'fig1','lick_psth');   % results directory
    lick_psth_summary_VP(vpcells,resdir_lick_psth,true)  % average lick PSTH
end

% Statistics on animal performance - FIG.1H
if lick_stat
    resdir_astat = fullfile(resdir,'fig1','anticipatory_stat');
    T13 = anticipatory_stat_VP('HDB13',resdir_astat);
    T17 = anticipatory_stat_VP('HDB17',resdir_astat);
    T25 = anticipatory_stat_VP('HDB25',resdir_astat);
    T32 = anticipatory_stat_VP('HDB32',resdir_astat);
    T38 = anticipatory_stat_VP('HDB38',resdir_astat);
    T = [{T13} {T17} {T25} {T32} {T38}];
    lineplot_with_errorbars_vp(T, {'HDB13' 'HDB17' 'HDB25' 'HDB32' 'HDB38'}, resdir_astat);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example PSTHs - FIG.2A-F
if response_example
    cellids = { 'HDB13_170228a_2.1';
        'HDB13_170301a_3.2';
        'HDB32_181121a_5.3';
        'HDB32_181130a_2.2';
        'HDB32_181130a_4.1';
        'HDB25_180401a_7.3';
        };
    resdir_raster_examples = fullfile(resdir,'fig2','raster_PSTH_examples');
    raster_PSTH_examples(cellids,resdir_raster_examples);
end

% Pie charts of responses - FIG.2G-I
if response_pie
    population = getcellresp(vpcells);   % output is used for Fig.2 pie charts
end

% Average PSTH and matrix plot - FIG.2J-L,O
if response_PSTH
    if ~response_pie
        population = getcellresp(vpcells);
    end
    resdir_avg_psth = fullfile(resdir,'fig2','avg_psth');
    avg_psth_VP(population.cue.excitation,population.cue.inhibition,'cueresponse',resdir_avg_psth); % FIG.2J-L - average response
    avg_psth_VP(population.reward.excitation,population.reward.inhibition,'rewardresponse',resdir_avg_psth);
    avg_psth_VP(population.punishment.excitation,population.punishment.inhibition,'punishresponse',resdir_avg_psth);
    
    % Plotting responses in a matrix - FIG.2O
    [~, G] = response_matrix(vpcells);
    saveas(G,fullfile(resdir_avg_psth,'imagesc_cellreponses.fig'))
    saveas(G,fullfile(resdir_avg_psth,'imagesc_cellreponses.eps'))
    
    % Compare response magnitude
    avg_psth_compare_maxval_VP(population.reward.excitation,'rewardresponse',population.punishment.excitation,'punishresponse','excitation',resdir_avg_psth);
    avg_psth_compare_maxval_VP(population.reward.inhibition,'rewardresponse',population.punishment.inhibition,'punishresponse','inhibition',resdir_avg_psth);
end

% Peak and latency statistics - Fig.2M-N
if box_stat
    latencies = latency_dist(vpcells,responsedir);  % latency properties of cue, reward and punishment response
    
    % Statistics, box plots
    resdir_boxstat = fullfile(resdir,'fig2','boxstat');
    if ~isfolder(resdir_boxstat)
        mkdir(resdir_boxstat)
    end
    [H1, Wp2] = boxstat(latencies.cue.excited.latency_peak,latencies.reward.excited.latency_peak,'cue','reward', 0.005);
    [H2, Wp2] = boxstat(latencies.cue.excited.latency_peak,latencies.punishment.excited.latency_peak,'cue','punishment', 0.005);
    [H3, Wp2] = boxstat(latencies.reward.excited.latency_peak,latencies.punishment.excited.latency_peak,'reward','punishment', 0.005);
    
    [H4, Wp2] = boxstat(latencies.cue.inhibited.latency_peak,latencies.reward.inhibited.latency_peak,'cue','reward', 0.005);
    [H5, Wp2] = boxstat(latencies.cue.inhibited.latency_peak,latencies.punishment.inhibited.latency_peak,'cue','punishment', 0.005);
    [H6, Wp2] = boxstat(latencies.reward.inhibited.latency_peak,latencies.punishment.inhibited.latency_peak,'reward','punishment', 0.005);
    
    % Save figures
    saveas(H1,fullfile(resdir_boxstat,'boxstat_excitation_cue_reward.fig'))
    set(H1,'Renderer','painters')
    saveas(H1,fullfile(resdir_boxstat,'boxstat_excitation_cue_reward.eps'))
    
    saveas(H2,fullfile(resdir_boxstat,'boxstat_excitation_cue_punishment.fig'))
    set(H2,'Renderer','painters')
    saveas(H2,fullfile(resdir_boxstat,'boxstat_excitation_cue_punishment.eps'))
    
    saveas(H3,fullfile(resdir_boxstat,'boxstat_excitation_reward_punishment.fig'))
    set(H3,'Renderer','painters')
    saveas(H3,fullfile(resdir_boxstat,'boxstat_excitation_reward_punishment.eps'))
    
    saveas(H4,fullfile(resdir_boxstat,'boxstat_inhibition_cue_reward.fig'))
    set(H4,'Renderer','painters')
    saveas(H4,fullfile(resdir_boxstat,'boxstat_inhibition_cue_reward.eps'))
    
    saveas(H5,fullfile(resdir_boxstat,'boxstat_inhibition_cue_punishment.fig'))
    set(H5,'Renderer','painters')
    saveas(H5,fullfile(resdir_boxstat,'boxstat_inhibition_cue_punishment.eps'))
    
    saveas(H6,fullfile(resdir_boxstat,'boxstat_inhibition_reward_punishment.fig'))
    set(H6,'Renderer','painters')
    saveas(H6,fullfile(resdir_boxstat,'boxstat_inhibition_reward_punishment.eps'))
    close all
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PSTH comparison of expected vs. unexpected - FIG.3A,C,E
if compare_exp_psth
    resdir_compare_exp_psth = fullfile(resdir,'fig3','compare_exp_psth');
    responsetype = [1, -1]; % excitation or onhibition
    for i = 1:length(responsetype)
        responsespec = responsetype(i);
        compare_expectations_VP(vpcells,'cueresponse',responsespec,resdir_compare_exp_psth)
        compare_expectations_VP(vpcells,'rewardresponse',responsespec,resdir_compare_exp_psth)
        compare_expectations_VP(vpcells,'punishresponse',responsespec,resdir_compare_exp_psth)
    end
end

% Boxplot comparison of expected vs. unexpected - FIG.3B,D,F
if compare_exp_box
    resdir_compare_exp_box = fullfile(resdir,'fig3','compare_exp_psth_box');
    responsetype = [1, -1]; %excitation or onhibition
    for i = 1:length(responsetype)
        responsespec = responsetype(i);
        compare_expectation_boxstat(vpcells,'cueresponse',responsespec,resdir_compare_exp_box)
        compare_expectation_boxstat(vpcells,'rewardresponse',responsespec,resdir_compare_exp_box)
        compare_expectation_boxstat(vpcells,'punishresponse',responsespec,resdir_compare_exp_box)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Autocorrelation examples - FIG.4A
if autocorr_bursting  % autocorrelation of neurons - to decide bursting, rhytmic etc..neurons
    segfilter = 'stim_excl_vp';
    filterinput = {'light_activation_duration',[-5 5],'margins',[0 0]};
    resdir_acg_burst = fullfile(resdir,'fig4','acg_burst');
    if ~isfolder(resdir_acg_burst)
        mkdir(resdir_acg_burst)
    end
    acg({'HDB17_170717a_1.4','HDB13_170306a_5.2',},0.5,'resdir',resdir_acg_burst,...
        'segfilter',segfilter,'filterinput',filterinput,...
        'minspikeno',100,'maxspikeno',10000,'issave',true); % Example ACGs
    if ~perform_acg
        acg_resdir = fullfile(resdir,'vpacg');   % results directory
    end
    resdir_sorting_burst = fullfile(resdir,'fig4','acg_sorted_burst');
    resdir_sorting_nonburst = fullfile(resdir,'fig4','acg_sorted_nonburst');
    acg_sorting(vpcells, acg_resdir, resdir_sorting_burst, resdir_sorting_nonburst) % sorting bursting and non-bursting neurons
end

if pie_PSTH_bursting
    
    % Pie charts - FIG.4B-C
    load(fullfile(getpref('cellbase','datapath'), '_paper_figs\code_review3\ACG_matrices.mat'));  % load ACG matrices (generated by vpacg.m)
    bursting_cells = cellids(BurstIndex>0.2);  % numbers for FIG.4B
    nonbursting_cells = cellids(BurstIndex<=0.2);
    
    bursting = getcellresp(bursting_cells); % numbers for FIG.4C
    nonbursting = getcellresp(nonbursting_cells);
    
    % PSTH averages - FIG.4D-I
    resdir_avg_psth_VP1 = fullfile(resdir,'fig4','bursting04');
    avg_psth_VP(bursting.cue.excitation,bursting.cue.inhibition,'cueresponse',resdir_avg_psth_VP1);
    avg_psth_VP(bursting.reward.excitation,bursting.reward.inhibition,'rewardresponse',resdir_avg_psth_VP1);
    avg_psth_VP(bursting.punishment.excitation,bursting.punishment.inhibition,'punishresponse',resdir_avg_psth_VP1);
    
    resdir_avg_psth_VP2 = fullfile(resdir,'fig4','nonbursting04');
    avg_psth_VP(nonbursting.cue.excitation,nonbursting.cue.inhibition,'cueresponse',resdir_avg_psth_VP2);
    avg_psth_VP(nonbursting.reward.excitation,nonbursting.reward.inhibition,'rewardresponse',resdir_avg_psth_VP2);
    avg_psth_VP(nonbursting.punishment.excitation,nonbursting.punishment.inhibition,'punishresponse',resdir_avg_psth_VP2);
    
    %      Compare responsive neuronal proportions - chi square-test
    [tbl1, chi2stat1, pval1] = chi2_vp(length([bursting.cue.excitation' bursting.cue.inhibition']),length(bursting_cells),length([nonbursting.cue.excitation' nonbursting.cue.inhibition']),length(nonbursting_cells));
    [tbl2, chi2stat2, pval2] = chi2_vp(length([bursting.reward.excitation' bursting.reward.inhibition']),length(bursting_cells),length([nonbursting.reward.excitation' nonbursting.reward.inhibition']),length(nonbursting_cells));
    [tbl3, chi2stat3, pval3] = chi2_vp(length([bursting.punishment.excitation' bursting.punishment.inhibition']),length(bursting_cells),length([nonbursting.punishment.excitation' nonbursting.punishment.inhibition']),length(nonbursting_cells));
    
    % Maxvalue stats
    resdir_maxval_comp_bursting = fullfile(resdir,'fig4','bursting_maxval04');
    avg_psth_compare_maxval_VP(bursting.cue.excitation,'cueresponse',nonbursting.cue.excitation,'cueresponse','excitation',resdir_maxval_comp_bursting);
    avg_psth_compare_maxval_VP(bursting.reward.excitation,'rewardresponse',nonbursting.reward.excitation,'rewardresponse','excitation',resdir_maxval_comp_bursting);
    avg_psth_compare_maxval_VP(bursting.punishment.excitation,'punishresponse',nonbursting.punishment.excitation, 'punishresponse','excitation',resdir_maxval_comp_bursting);
    
    avg_psth_compare_maxval_VP(bursting.cue.inhibition,'cueresponse',nonbursting.cue.inhibition,'cueresponse','inhibition',resdir_maxval_comp_bursting);
    avg_psth_compare_maxval_VP(bursting.reward.inhibition,'rewardresponse',nonbursting.reward.inhibition,'rewardresponse','inhibition',resdir_maxval_comp_bursting);
    avg_psth_compare_maxval_VP(bursting.punishment.inhibition,'punishresponse',nonbursting.punishment.inhibition,'punishresponse','inhibition',resdir_maxval_comp_bursting);
    
    % Average ACGs
    segfilter = 'stim_excl_vp'; % exclude recording segments with stimulation
    filterinput = {'light_activation_duration',[-5 5],'margins',[0 0]};
    resdir_acg_average = fullfile(resdir,'fig4','acg_average');
    if ~isfolder(resdir_acg_average)
        mkdir(resdir_acg_average)
    end
    avg_vpacg(bursting_cells,0.5,'resdir',resdir_acg_average,...
        'segfilter',segfilter,'filterinput',filterinput,...
        'minspikeno',100,'maxspikeno',10000,'issave',true);
    avg_vpacg(nonbursting_cells,0.5,'resdir',resdir_acg_average,...
        'segfilter',segfilter,'filterinput',filterinput,...
        'minspikeno',100,'maxspikeno',10000,'issave',true);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Autocorrelation examples - FIG.5A
if autocorr_rhythmic  % autocorrelation of neurons - to decide bursting, rhytmic etc..neurons
    segfilter = 'stim_excl_vp';
    filterinput = {'light_activation_duration',[-5 5],'margins',[0 0]};
    resdir_acg_example = fullfile(resdir,'fig5','acg_example');
    if ~isfolder(resdir_acg_example)
        mkdir(resdir_acg_example)
    end
    acg({'HDB25_180426a_8.2'},0.5,'resdir',resdir_acg_example,...
        'segfilter',segfilter,'filterinput',filterinput,...
        'minspikeno',100,'maxspikeno',10000,'issave',true); % Example ACGs
end

if pie_PSTH_rhythmic
    
    % Pie charts - FIG.5B-D
    resdir_acg_rhythmic = fullfile(resdir,'fig5','acg_rhythmic');  % save acg
    resdir2 = fullfile(resdir,'fig5','acg_nonrhythmic');  % save acg
    resdir3 = fullfile(resdir,'fig5','psth_rhythmic');  % for testing
    
    % Rhytmicity - calculate Beta and Gamma Index
    [phasic_cells, BetaIndex, GammaIndex, isbeta, isgamma] = pr_freq_detection_VP(vpcells,resdir_acg_rhythmic,resdir2,resdir3,acgdir,responsedir);
    beta_cells = phasic_cells(isbeta==1); % beta cells
    gamma_cells = phasic_cells(isgamma==1); % gamma cells
    
    rhythmic_cells = [beta_cells gamma_cells];  % rhythmic neurons
    rhythmic_cells = unique(rhythmic_cells);   % beta/gamma may overlap
    nonrhythmic_cells = setdiff(vpcells,rhythmic_cells);  % non-rhythmic VP neurons
    
    phasic = getcellresp(rhythmic_cells);  % FIG.5B-D
    nonphasic = getcellresp(nonrhythmic_cells);
    
    % PSTH averages - FIG.5E-J
    resdir_avg_psth_VP1 = fullfile(resdir,'fig5','phasic_avg_psth');
    avg_psth_VP(phasic.cue.excitation,phasic.cue.inhibition,'cueresponse',resdir_avg_psth_VP1);
    avg_psth_VP(phasic.reward.excitation,phasic.reward.inhibition,'rewardresponse',resdir_avg_psth_VP1);
    avg_psth_VP(phasic.punishment.excitation,phasic.punishment.inhibition,'punishresponse',resdir_avg_psth_VP1);
    
    resdir_avg_psth_VP2 = fullfile(resdir,'fig5','nonphasic_avg_psth');
    avg_psth_VP(nonphasic.cue.excitation,nonphasic.cue.inhibition,'cueresponse',resdir_avg_psth_VP2);
    avg_psth_VP(nonphasic.reward.excitation,nonphasic.reward.inhibition,'rewardresponse',resdir_avg_psth_VP2);
    avg_psth_VP(nonphasic.punishment.excitation,nonphasic.punishment.inhibition,'punishresponse',resdir_avg_psth_VP2);
    
    % Compare proportions of synchronous neurons in the two groups
    [tbl5,chi2stat5,pval5] = chi2_vp(length([phasic.cue.excitation' phasic.cue.inhibition']),length(rhythmic_cells),length([nonphasic.cue.excitation' nonphasic.cue.inhibition']),length(nonrhythmic_cells));
    [tbl6,chi2stat6,pval6] = chi2_vp(length([phasic.reward.excitation' phasic.reward.inhibition']),length(rhythmic_cells),length([nonphasic.reward.excitation' nonphasic.reward.inhibition']),length(nonrhythmic_cells));
    [tbl7,chi2stat7,pval7] = chi2_vp(length([phasic.punishment.excitation' phasic.punishment.inhibition']),length(rhythmic_cells),length([nonphasic.punishment.excitation' nonphasic.punishment.inhibition']),length(nonrhythmic_cells));
    
    % Maxvalue statistics
    resdir_maxval_comp_phasic = fullfile(resdir,'fig5','phasic_maxval');
    avg_psth_compare_maxval_VP(phasic.cue.excitation,'cueresponse',nonphasic.cue.excitation,'cueresponse','excitation',resdir_maxval_comp_phasic);
    avg_psth_compare_maxval_VP(phasic.reward.excitation,'rewardresponse',nonphasic.reward.excitation,'rewardresponse','excitation',resdir_maxval_comp_phasic);
    avg_psth_compare_maxval_VP(phasic.punishment.excitation,'punishresponse',nonphasic.punishment.excitation,'punishresponse','excitation',resdir_maxval_comp_phasic);
    
    avg_psth_compare_maxval_VP(phasic.cue.inhibition, 'cueresponse', nonphasic.cue.inhibition,'cueresponse','inhibition',resdir_maxval_comp_phasic);
    avg_psth_compare_maxval_VP(phasic.reward.inhibition, 'rewardresponse', nonphasic.reward.inhibition,'rewardresponse','inhibition',resdir_maxval_comp_phasic);
    avg_psth_compare_maxval_VP(phasic.punishment.inhibition, 'punishresponse', nonphasic.punishment.inhibition,'punishresponse','inhibition',resdir_maxval_comp_phasic);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Raster and PSTH for single spikes - FIG.6A,C,D,F
if single_spikes
    burst_detection(vpcells, 'single') % detecting single spikes
    cells = {'HDB13_170228a_2.1' , 'HDB13_170301a_3.1'};  % example neurons
    for i = 1:length(cells) % example PSTHs
        viewcell2p(cells{i},'TriggerName','StimulusOn','SortEvent','TrialStart','sigma',0.07,... % single spikes
            'eventtype','behav','ShowEvents',{{'DeliverAllFeedback'}},'Partitions','#TrialType' , 'window', [-5 5]);
        viewcell2b(cells{i},'TriggerName','StimulusOn','SortEvent','TrialStart','sigma',0.07,... % all spikes
            'eventtype','behav','ShowEvents',{{'DeliverAllFeedback'}},'Partitions','#TrialType','window',[-5 5]);
    end
end

% Raster and PSTH for bursts - FIG.6B,E
if bursts
    burst_detection(vpcells, 'burst') % detecting bursts
    cells = {'HDB13_170228a_2.1' , 'HDB13_170301a_3.1'};   % example neurons
    for i = 1:length(cells)  % example PSTHs
        viewcell2b(cells{i},'TriggerName','StimulusOn','SortEvent','TrialStart','sigma',0.07,... % all spikes
            'eventtype','behav','ShowEvents',{{'DeliverAllFeedback'}},'Partitions','#TrialType','window',[-5 5]);
        viewcell2burst(cells{i},'TriggerName','StimulusOn','SortEvent','TrialStart','sigma',0.07,... % burst spikes
            'eventtype','behav','ShowEvents',{{'DeliverAllFeedback'}},'Partitions','#TrialType','window',[-5 5]);
    end
end

% Find cells with differential burst and single spike response - FIG.6G-H
if pie_burst_single 
    resdir_burst_vs_spike = fullfile(resdir,'fig6','burst_vs_single_avg_psth'); % results directory
    resdir_pie6 = fullfile(resdir,'fig6');
    pval = 0.01;
    cond = {'cueresponse' 'rewardresponse' 'punishresponse'};
    for k = 1:length(cond)
        [besiinx_cell, biseinx_cell] = redo_vpresponsesorter(vpcells, cond{k}, pval, resdir_pie6); % find those cells
        avg_psth_VP(besiinx_cell,besiinx_cell,cond{k},fullfile(resdir_burst_vs_spike,'excitation_b'),'burst'); % average PSTH (burst activation, single spike inhibition)
        avg_psth_VP(biseinx_cell,biseinx_cell,cond{k},fullfile(resdir_burst_vs_spike,'excitation_s'),'burst'); % average PSTH (single spike activation, burst inhibition) - FIG.S6
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Monosynaptic pairs - FIG.7A-B
if select_pairs
    most_pairs = select_monosyn_pairs(ccgdir); % only part of the connections are shown in the fig.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if vpcross
    
    % Example cross_correlogram - FIG.8A
    resdir_nontetrode = fullfile(resdir,'fig8','vpccg','nontetrodepairs');
    if ~isfolder(resdir_nontetrode)
        mkdir(resdir_nontetrode)
    end
    segfilter = 'stim_excl_vp';
    filterinput = {'light_activation_duration',[-5 5],'margins',[0 0]};
    ccg({'HDB17_170718a_2.2','HDB17_170718a_4.1'},0.05,'whichcells','nontetrodepairs','resdir',resdir_nontetrode,...
        'segfilter',segfilter,'filterinput',filterinput,...
        'minspikeno',100,'maxspikeno',10000,'issave',true);
    
    % Average cross correlogram - FIG.8B
    avg_vpccg;
    saveas(gcf,fullfile(ccgdir,'avg_ccg.fig'))
    set(gcf,'Renderer','painters')
    saveas(gcf,fullfile(ccgdir,'avg_ccg.eps'))
end

if sync_pie_PSTH
    
    % Pie chart - FIG.8C
    [sync_cellids, nonsync_cellids] = sync_pair_responses(vpcells,ccgdir); % synchronous and nonsynchronous cell groups (for pie chart)
    sync = getcellresp(sync_cellids); % synchronous cellgroups
    nonsync = getcellresp(nonsync_cellids); % nonsynchronous cellgroups
    
    % Average PSTH - sync - FIG.8D-F
    resdir_sync = fullfile(resdir,'fig8','sync_avg_PSTH');
    avg_psth_VP(sync.cue.excitation,sync.cue.inhibition,'cueresponse',resdir_sync);
    avg_psth_VP(sync.reward.excitation,sync.reward.inhibition,'rewardresponse',resdir_sync);
    avg_psth_VP(sync.punishment.excitation,sync.punishment.inhibition,'punishresponse',resdir_sync);
    
    % Average PSTH - nonsync - FIG.8G-I
    resdir_nonsync = fullfile(resdir,'fig8','nonsync_avg_PSTH');
    avg_psth_VP(nonsync.cue.excitation,nonsync.cue.inhibition,'cueresponse',resdir_nonsync);
    avg_psth_VP(nonsync.reward.excitation,nonsync.reward.inhibition,'rewardresponse',resdir_nonsync);
    avg_psth_VP(nonsync.punishment.excitation,nonsync.punishment.inhibition,'punishresponse',resdir_nonsync);
    
    % Compare proportions of synchronous neurons in the two groups
    [tbl8,chi2stat8,pval8] = chi2_vp(length([sync.cue.excitation' sync.cue.inhibition']), length(sync_cellids),length([nonsync.cue.excitation' nonsync.cue.inhibition']), length(nonsync_cellids));
    [tbl9,chi2stat9,pval9] = chi2_vp(length([sync.reward.excitation' sync.reward.inhibition']), length(sync_cellids),length([nonsync.reward.excitation' nonsync.reward.inhibition']), length(nonsync_cellids));
    [tbl10,chi2stat10,pval10] = chi2_vp(length([sync.punishment.excitation' sync.punishment.inhibition']), length(sync_cellids),length([nonsync.punishment.excitation' nonsync.punishment.inhibition']), length(nonsync_cellids));
    
    % Maxvalue statistics
    resdir_maxval_comp_sync = fullfile(resdir,'fig8','sync_maxval');
    avg_psth_compare_maxval_VP(sync.cue.excitation,'cueresponse',nonsync.cue.excitation,'cueresponse','excitation',resdir_maxval_comp_sync);
    avg_psth_compare_maxval_VP(sync.reward.excitation,'rewardresponse',nonsync.reward.excitation,'rewardresponse','excitation',resdir_maxval_comp_sync);
    avg_psth_compare_maxval_VP(sync.punishment.excitation,'punishresponse',nonsync.punishment.excitation,'punishresponse','excitation',resdir_maxval_comp_sync);
    
    avg_psth_compare_maxval_VP(sync.cue.inhibition,'cueresponse',nonsync.cue.inhibition,'cueresponse','inhibition',resdir_maxval_comp_sync);
    avg_psth_compare_maxval_VP(sync.reward.inhibition,'rewardresponse',nonsync.reward.inhibition,'rewardresponse','inhibition',resdir_maxval_comp_sync);
    avg_psth_compare_maxval_VP(sync.punishment.inhibition,'punishresponse',nonsync.punishment.inhibition,'punishresponse','inhibition',resdir_maxval_comp_sync);
    
end

if norhythmandburst  % FIG.8K - select non-bursting- non-rhythmic cells and compare their synchrony to others
    load(fullfile(acgdir,'ACG_matrices.mat'));  % load ACG matrices
    nonbursting_cells = cellids(BurstIndex<=0.2); % select nonbursting cells
    
    if ~pie_PSTH_rhythmic
        resdir_acg_rhythmic = fullfile(resdir,'fig5','acg_rhythmic');  % rhythmic acg folder
    end
    % Rhytmicity - calculate Beta and Gamma Index
    resdir_acg_rhythmic = fullfile(resdir,'fig5','acg_rhythmic');  % save acg
    resdir2 = fullfile(resdir,'fig5','acg_nonrhythmic');  % save acg
    resdir3 = fullfile(resdir,'fig5','psth_rhythmic');  % for testing
    [phasic_cells, BetaIndex, GammaIndex, isbeta, isgamma] = pr_freq_detection_VP(vpcells,resdir_acg_rhythmic,resdir2,resdir3,acgdir,responsedir);
    beta_cells = phasic_cells(isbeta==1); % beta cells
    gamma_cells = phasic_cells(isgamma==1); % gamma cells
    
    rhythmic_cells = [beta_cells gamma_cells];  % rhythmic neurons
    rhythmic_cells = unique(rhythmic_cells);   % beta/gamma may overlap
    nonrhythmic_cells = setdiff(vpcells,rhythmic_cells);  % non-rhythmic VP neurons
    nonryhtmic_nonburst = intersect(nonbursting_cells, nonrhythmic_cells); % non-rhythmic-nonbursting cells
    rest = setdiff(vpcells, nonryhtmic_nonburst);
    
    % Sorting non-rhythmic-nonbursting neurons and the rest according to synchrony
    [sync_cellids1, nonsync_cellids1] = sync_pair_responses(nonryhtmic_nonburst,ccgdir); % synchronous and nonsynchronous cell groups (for pie chart)
    [sync_cellids2, nonsync_cellids2] = sync_pair_responses(rest,ccgdir); % synchronous and nonsynchronous cell groups (for pie chart)
    
    % Compare proportions of synchronous neurons in the two groups
    [tbl4,chi2stat4,pval4] = chi2_vp(length(sync_cellids1), length(nonryhtmic_nonburst),length(sync_cellids2), length(rest));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUPPLEMENTARY ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

animals = {'HDB13' 'HDB17' 'HDB25' 'HDB32' 'HDB38'}; % animals involved in analysis
if figS1spec  % Fig.S1
    
    % Neurons per animal
    [all_cells, cells_overall] = peranimal(vpcells,animals);
end

if figS2spec  % Fig.S2
    if ~figS1spec % if neurons are not grouped by animals
        
        % Neurons per animal
        [all_cells, cells_overall] = peranimal(vpcells,animals);
    end
    
    for k = 1:length(animals) % loop through animals - FIGS2
        resdir_compare_exp_psth_S = fullfile(resdir,'figS2','compare_exp_psth', animals{k});
        resdir_compare_exp_box_S = fullfile(resdir,'figS2','compare_exp_psth_box', animals{k});
        responsetype = [-1, 1]; %excitation or onhibition
        for i = 1:length(responsetype)
            responsespec = responsetype(i);
            compare_expectations_VP(all_cells.(animals{k}),'cueresponse',responsespec,resdir_compare_exp_psth_S) % average PSTH
            compare_expectation_boxstat(all_cells.(animals{k}),'cueresponse',responsespec,resdir_compare_exp_box_S) % boxstat
        end
    end
end

if figS3spec  % Fig.S3
    resdir_s3 = fullfile(resdir,'figS3');
    [bursting_per_animal, nonbursting_per_animal, b_cellids_per_animal, nonb_cellids_per_animal] = cell_type_proportions(vpcells, animals, 'bursting', resdir_s3, resdir);
    [b_response_per_animal] = sortresponses_animal(animals, b_cellids_per_animal);
    [nb_response_per_animal] = sortresponses_animal(animals, nonb_cellids_per_animal);
end

if figS4spec  % Fig.S4
    resdir_s4 = fullfile(resdir,'figS4');
    [rhythmic_per_animal, nonrhythmic_per_animal, r_cellids_per_animal, nonr_cellids_per_animal] = cell_type_proportions(vpcells, animals, 'rhythmic', resdir_s4, resdir);
    [r_response_per_animal] = sortresponses_animal(animals, r_cellids_per_animal);
    [nr_response_per_animal] = sortresponses_animal(animals, nonr_cellids_per_animal);
end

if figS5spec  % Fig.S5
    resdir_s5 = fullfile(resdir,'figS5');
    [synchronous_per_animal, asynchronous_per_animal, s_cellids_per_animal, as_cellids_per_animal] = cell_type_proportions(vpcells, animals, 'sync', resdir_s5, resdir);
    [s_response_per_animal] = sortresponses_animal(animals, s_cellids_per_animal);
    [as_response_per_animal] = sortresponses_animal(animals, as_cellids_per_animal);
end
end

% -------------------------------------------------------------------------
function responses = getcellresp(cellids)

% VP neuron responses
cueresp = getvalue('cueresponse',cellids);
rewardresp = getvalue('rewardresponse',cellids);
punishmentresp = getvalue('punishresponse',cellids);

% Cue response
cue_e = cellids(cueresp == 1);  % activated
cue_i = cellids(cueresp == -1); % inhibited
cue_n = cellids(cueresp == 0);  % non-responsive

%reward response
rew_e = cellids(rewardresp == 1);  % activated
rew_i = cellids(rewardresp == -1); % inhibited
rew_n = cellids(rewardresp == 0);  % non-responsive

%punishment response
pun_e = cellids(punishmentresp == 1);  % activated
pun_i = cellids(punishmentresp == -1); % inhibited
pun_n = cellids(punishmentresp == 0);  % non-responsive

% Output
responses = struct;
responses.cue.excitation = cue_e';
responses.cue.inhibition = cue_i';
responses.cue.none = cue_n';
responses.reward.excitation = rew_e';
responses.reward.inhibition = rew_i';
responses.reward.none = rew_n';
responses.punishment.excitation = pun_e';
responses.punishment.inhibition = pun_i';
responses.punishment.none = pun_n';
end

% -------------------------------------------------------------------------
function [all_cells, cells_overall] = peranimal(vpcells,animals)

% Neurons per animal
cells_overall = struct;
all_cells = struct;
for g = 1:length(animals)
    aID = animals{g};
    cells_per_animal = [];
    for h = 1:length(vpcells)
        [ratname,~,~,~] = cellid2tags(vpcells{h});
        if strcmp(ratname,aID)
            cells_per_animal = [cells_per_animal {vpcells{h}}]; % append cell
        end
    end
    all_cells.(aID) = cells_per_animal;
    cells_overall.(aID) = getcellresp(cells_per_animal); % get cell responses for eact animal (FIGS1.C)
end
end