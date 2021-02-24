function slowanalysis_pavlovian(vpcells, preprocess, fig1spec, fig2spec, ...
    fig3spec, fig4spec, fig5spec, fig6spec, fig7spec, fig8spec, suppl_figspec, rev_figspec)
%SLOWANALYSIS_PAVLOVIAN   Analysis main function for VP project.
%   SLOWANALYSIS_PAVLOVIAN(VPCELLS, PRE, F1, F2, ..., F8, S1, R1) is the main
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
%   Code for review process (iScience) added 18/02/2021

% Input argument check
narginchk(0,12)
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
vpcross = fig8spec(1); % FIG.8A-B - example and average CCG
sync_pie_PSTH = fig8spec(2); % FIG.8C-I - pie chart and average PSTH for sync and non-sync
norhythmandburst = fig8spec(3); % FIG.8K - pie chart of sync firing in the non-bursting-non-rhythmic group

if nargin < 11  % do supplementary analysis or not
    suppl_figspec = [1 1 1 1 1];
end
figS1spec = suppl_figspec(1); % FigS1
figS2spec = suppl_figspec(2); % FigS2
figS3spec = suppl_figspec(3); % FigS3
figS4spec = suppl_figspec(4); % FigS4
figS5spec = suppl_figspec(5); % FigS5

if nargin < 12  % figures for revision
    rev_figspec = [1 1 1 1 1 1 1 1 1];
end
figR1spec = rev_figspec(1);
figR2spec = rev_figspec(2);
figR3spec = rev_figspec(3);
figR4spec = rev_figspec(4);
figR5spec = rev_figspec(5);
figR6spec = rev_figspec(6);
figR8spec = rev_figspec(7);
figR9spec = rev_figspec(8);
figR10spec = rev_figspec(9);

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

if figR6spec % Exclude possible replicates of recorded neurons- FIGR6
    Cells2Exclude = trackoverdays(vpcells);
    vpcells = setdiff(vpcells, Cells2Exclude);
end

if figR9spec % FIGR9 - exclude PV-Cre animal
    VP17cells = findcell('Rat', 'HDB17')
    vpcells = setdiff(vpcells, VP17cells);
end

% Grouping neurons based on their response during behavior - properties added to TheMatrix
response_resdir = fullfile(resdir,'vpresponsesorter');   % results directory
if response_analysis
    vpresponsesorter(vpcells,1,response_resdir,'cue');
    vpresponsesorter(vpcells,1,response_resdir,'rew');
    vpresponsesorter(vpcells,1,response_resdir,'pun');
    vpresponsesorter(vpcells,1,response_resdir,'om');
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
    resdir_compare_exp_psth = fullfile(resdir,'fig3_allcells','compare_exp_psth');
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
    resdir_avg_psth_VP1 = fullfile(resdir,'fig4','bursting30');
    avg_psth_VP(bursting.cue.excitation,bursting.cue.inhibition,'cueresponse',resdir_avg_psth_VP1);
    avg_psth_VP(bursting.reward.excitation,bursting.reward.inhibition,'rewardresponse',resdir_avg_psth_VP1);
    avg_psth_VP(bursting.punishment.excitation,bursting.punishment.inhibition,'punishresponse',resdir_avg_psth_VP1);
    
    resdir_avg_psth_VP2 = fullfile(resdir,'fig4','nonbursting30');
    avg_psth_VP(nonbursting.cue.excitation,nonbursting.cue.inhibition,'cueresponse',resdir_avg_psth_VP2);
    avg_psth_VP(nonbursting.reward.excitation,nonbursting.reward.inhibition,'rewardresponse',resdir_avg_psth_VP2);
    avg_psth_VP(nonbursting.punishment.excitation,nonbursting.punishment.inhibition,'punishresponse',resdir_avg_psth_VP2);
    
    % Compare responsive neuronal proportions - chi square-test
    [tbl1, chi2stat1, pval1] = chi2_vp(length([bursting.cue.excitation' bursting.cue.inhibition']),length(bursting_cells),length([nonbursting.cue.excitation' nonbursting.cue.inhibition']),length(nonbursting_cells));
    [tbl2, chi2stat2, pval2] = chi2_vp(length([bursting.reward.excitation' bursting.reward.inhibition']),length(bursting_cells),length([nonbursting.reward.excitation' nonbursting.reward.inhibition']),length(nonbursting_cells));
    [tbl3, chi2stat3, pval3] = chi2_vp(length([bursting.punishment.excitation' bursting.punishment.inhibition']),length(bursting_cells),length([nonbursting.punishment.excitation' nonbursting.punishment.inhibition']),length(nonbursting_cells));
    
    % Maxvalue stats
    resdir_maxval_comp_bursting = fullfile(resdir,'fig4','bursting_maxval30');
    avg_psth_compare_maxval_VP(bursting.cue.excitation,'cueresponse',nonbursting.cue.excitation,'cueresponse','excitation',resdir_maxval_comp_bursting);
    avg_psth_compare_maxval_VP(bursting.reward.excitation,'rewardresponse',nonbursting.reward.excitation,'rewardresponse','excitation',resdir_maxval_comp_bursting);
    avg_psth_compare_maxval_VP(bursting.punishment.excitation,'punishresponse',nonbursting.punishment.excitation, 'punishresponse','excitation',resdir_maxval_comp_bursting);
    
    avg_psth_compare_maxval_VP(bursting.cue.inhibition,'cueresponse',nonbursting.cue.inhibition,'cueresponse','inhibition',resdir_maxval_comp_bursting);
    avg_psth_compare_maxval_VP(bursting.reward.inhibition,'rewardresponse',nonbursting.reward.inhibition,'rewardresponse','inhibition',resdir_maxval_comp_bursting);
    avg_psth_compare_maxval_VP(bursting.punishment.inhibition,'punishresponse',nonbursting.punishment.inhibition,'punishresponse','inhibition',resdir_maxval_comp_bursting);
    
    % Average ACGs
    segfilter = 'stim_excl_vp'; % exclude recording segments with stimulation
    filterinput = {'light_activation_duration',[-5 5],'margins',[0 0]};
    resdir_acg_average = fullfile(resdir,'fig4','acg_average30');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Analysis for review @ iScience %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bursting and non-bursting cells
load(fullfile(getpref('cellbase','datapath'), '_paper_figs\code_review3\ACG_matrices.mat'));  % load ACG matrices (generated by vpacg.m)
bursting_cells = cellids(BurstIndex>0.2);  % numbers for FIG.4B
nonbursting_cells = cellids(BurstIndex<=0.2);

% Rhythmic and non-rhythmic cells
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

% Synchronous and asynchronous cells
[sync_cellids, nonsync_cellids] = sync_pair_responses(vpcells,ccgdir); % synchronous and nonsynchronous cell groups (for pie chart)
sync = getcellresp(sync_cellids); % synchronous cellgroups
nonsync = getcellresp(nonsync_cellids); % nonsynchronous cellgroups

rootdir = fullfile(getpref('cellbase', 'datapath'), '_paper_figs\Revision_iScience'); % all additional analysis is saved here

if figR1spec
    
    % Show that airpuff is way louder than just the valve click -FIGR1.B
    sourcedir = 'F:\auditory_pavlovian_cellbase\_paper_figs\Revision_iScience\_analysis\';
    [valve_data, header, ~] = xlsread(fullfile(sourcedir, 'valve_vs_airpuff.xlsx'));
    median_airpuff = median(valve_data(:,1));
    median_valve = median(valve_data(:,2));
    se_airpuff = se_of_median(valve_data(:,1));
    se_valve = se_of_median(valve_data(:,2));
    [h, p] = boxstat(valve_data(:,1),valve_data(:,2), 'airpuff', 'valve', 0.05)
    figure;
    bar([1 2], [median_airpuff, median_valve])
    hold on
    errorbar([1 2],  [median_airpuff, median_valve], [-se_airpuff, -se_valve], [se_airpuff, se_valve]);
    set(gcf, 'renderer', 'painters')
    saveas(gcf,fullfile(resdir_R1, 'airpuff_vs_valve.eps'));
    
    % Show the lack of neuronal response when omission occurs % - FIGR1.C
    resdir_R1 = fullfile(rootdir, '_analysis', 'R1'); % results directories
    resdir_avg_psth_o = fullfile(resdir_R1,'omissionresponse'); % results directory for omission 'responsive' cells aligned to feedback
    resdir_avg_psth_o = fullfile(resdir_R1,'omissionresponse_rewardcells'); % results directory for reward responsive cells aligned to feedback
    resdir_avg_psth_o = fullfile(resdir_R1,'omissionresponse_punishcells'); % results directory for punishment responsive cells aligned to feedback
    
    population = getcellresp(vpcells);   % output is used for Fig.R1 pie charts
    avg_psth_VP(population.omission.excitation,population.omission.inhibition,'omissionresponse',resdir_avg_psth_o); % average PSTH
    
    % Reward and punishment responsive cells are not responding to the lack of reinforcement
    avg_psth_VP(population.reward.excitation,population.reward.inhibition,'omissionresponse',resdir_avg_psth_o); % FigR1.D
    avg_psth_VP(population.punishment.excitation,population.punishment.inhibition,'omissionresponse',resdir_avg_psth_o);
end

if figR2spec % analize bursting neurons using a 30ms ISI cutoff
    acg_resdir30 = fullfile(rootdir,'vpacg30');   % results directory
    vpacg30(vpcells,acg_resdir30,true); % autocorrelograms and burst indices are recalculated
    
    % Data used for pie charts - FIGR2.A-B
    load(fullfile(acg_resdir30, 'ACG_matrices.mat'));  % load ACG matrices (generated by vpacg30.m)
    bursting_cells30 = cellids(BurstIndex>0.2); % cellids of bursting cells
    nonbursting_cells30 = cellids(BurstIndex<=0.2); % cellids of non-bursting cells
    
    bursting30 = getcellresp(bursting_cells30); % cue, reward and punishment response of bursting cells
    nonbursting30 = getcellresp(nonbursting_cells30); % cue, reward and punishment response of non-bursting cells
    
    % PSTH averages - FIGR2.C-H
    resdir_avg_psth_R2_1 = fullfile(rootdir, '_analysis', 'R2', 'bursting30');
    avg_psth_VP(bursting30.cue.excitation,bursting30.cue.inhibition,'cueresponse',resdir_avg_psth_R2_1); %C
    avg_psth_VP(bursting30.reward.excitation,bursting30.reward.inhibition,'rewardresponse',resdir_avg_psth_R2_1); %D
    avg_psth_VP(bursting30.punishment.excitation,bursting30.punishment.inhibition,'punishresponse',resdir_avg_psth_R2_1); %E
    
    resdir_avg_psth_R2_2 = fullfile(rootdir, '_analysis', 'R2', 'nonbursting30');
    avg_psth_VP(nonbursting30.cue.excitation,nonbursting30.cue.inhibition,'cueresponse',resdir_avg_psth_R2_2); %F
    avg_psth_VP(nonbursting30.reward.excitation,nonbursting30.reward.inhibition,'rewardresponse',resdir_avg_psth_R2_2); %G
    avg_psth_VP(nonbursting30.punishment.excitation,nonbursting30.punishment.inhibition,'punishresponse',resdir_avg_psth_R2_2); %H
    
    % Compare responsive neuronal proportions across bursting and non-bursting cells - chi square-test
    [tblR2_1, chi2statR2_1, pvalR2_1] = chi2_vp(length([bursting30.cue.excitation' bursting30.cue.inhibition']),length(bursting_cells30),length([nonbursting30.cue.excitation' nonbursting30.cue.inhibition']),length(nonbursting_cells30));
    [tblR2_2, chi2statR2_2, pvalR2_2] = chi2_vp(length([bursting30.reward.excitation' bursting30.reward.inhibition']),length(bursting_cells30),length([nonbursting30.reward.excitation' nonbursting30.reward.inhibition']),length(nonbursting_cells30));
    [tblR2_3, chi2statR2_3, pvalR2_3] = chi2_vp(length([bursting30.punishment.excitation' bursting30.punishment.inhibition']),length(bursting_cells30),length([nonbursting30.punishment.excitation' nonbursting30.punishment.inhibition']),length(nonbursting_cells30));
    
    % Maxvalue statistics
    resdir_avg_psth_R2_3 = fullfile(rootdir, '_analysis', 'R2', 'bursting_maxval30');
    avg_psth_compare_maxval_VP(bursting30.cue.excitation,'cueresponse',nonbursting30.cue.excitation,'cueresponse','excitation',resdir_avg_psth_R2_3);
    avg_psth_compare_maxval_VP(bursting30.reward.excitation,'rewardresponse',nonbursting30.reward.excitation,'rewardresponse','excitation',resdir_avg_psth_R2_3);
    avg_psth_compare_maxval_VP(bursting30.punishment.excitation,'punishresponse',nonbursting30.punishment.excitation, 'punishresponse','excitation',resdir_avg_psth_R2_3);
    
    avg_psth_compare_maxval_VP(bursting30.cue.inhibition,'cueresponse',nonbursting30.cue.inhibition,'cueresponse','inhibition',resdir_avg_psth_R2_3);
    avg_psth_compare_maxval_VP(bursting30.reward.inhibition,'rewardresponse',nonbursting30.reward.inhibition,'rewardresponse','inhibition',resdir_avg_psth_R2_3);
    avg_psth_compare_maxval_VP(bursting30.punishment.inhibition,'punishresponse',nonbursting30.punishment.inhibition,'punishresponse','inhibition',resdir_avg_psth_R2_3);
end

if figR3spec % Spike shape analysis
    if ~ismember('SpikeShape',listtag('prop')) % perform spike shape analysis if properties are not added to CellBase yet
        resdirR3_1 = fullfile(getpref('cellbase', 'datapath'),'spikeshapeanalysis_newdata'); % saveresults here
        vp_spikeshape_analysis(resdirR3_1);
    end
    
    % Compare spike shape properties of electrophysiological groups - FIGR3 boxplots
    compare_spikeshape([{bursting_cells} {nonbursting_cells}], [{'bursting'} {'nonbursting'}] ,fullfile(rootdir,'_analysis', 'R3', 'bursting_spikeanalysis'));
    compare_spikeshape([{rhythmic_cells} {nonrhythmic_cells}], [{'rhythmic'}, {'nonrhythmic'}] ,fullfile(rootdir,'_analysis', 'R3', 'rhythmic_spikeanalysis'));
    compare_spikeshape([{sync_cellids} {nonsync_cellids}], [{'synchronous'}, {'asynchronous'}] ,fullfile(rootdir,'_analysis', 'R3', 'sync_spikeanalysis'));
    
    resdirR3_2 = fullfile(rootdir, '_analysis', 'R3', 'avg_spikeshape'); % results directory
    avg_spikeshape(bursting_cells, resdirR3_2, 'bursting'); % average waveform of bursting neurons
    avg_spikeshape(nonbursting_cells, resdirR3_2, 'nonbursting'); % average waveform of non-bursting neurons
    
    avg_spikeshape(rhythmic_cells, resdirR3_2, 'rhythmic'); % average waveform of rhythmic neurons
    avg_spikeshape(nonrhythmic_cells, resdirR3_2, 'nonrhythmic'); % average waveform of non-rhythmic neurons
    
    avg_spikeshape(sync_cellids, resdirR3_2, 'synchronous'); % average waveform of synchronous neurons
    avg_spikeshape(nonsync_cellids, resdirR3_2, 'asynchronous'); % average waveform of asynchronous neurons
end

if figR4spec % Connection between bursting and firing rate - FIGR4.A-B
    if ~ismember('baseline_FR',listtag('prop'))
        resdirR4_1 = fullfile(rootdir, '_analysis', 'R4', 'baseline_FR');
        extract_baselineFR(vpcells,'cue',resdirR4_1); % Calculate baseline firing rate and add it to CellBase
    end
    resdirR4_2 = fullfile(rootdir, '_analysis', 'R4', 'baseline_vs_burst');
    baseline_vs_burst(vpcells, resdirR4_2); %FIGR4A-B
end

if figR5spec % Firing rate distribution of VP e-types
    if ~ismember('baseline_FR',listtag('prop'))
        resdirR4_1 = fullfile(rootdir, '_analysis', 'R4', 'baseline_FR');
        extract_baselineFR(vpcells,'cue',resdirR4_1); % Calculate baseline firing rate and add it to CellBase
    end
    
    BFR = getvalue('baseline_FR', vpcells); % Get baseline firing rate from The Matrix
    
    [~, b_cellinx, ~] = intersect(vpcells, bursting_cells); %FIGR5.A
    [~, nb_cellinx, ~] = intersect(vpcells, nonbursting_cells);
    boxstat(BFR(b_cellinx), BFR(nb_cellinx), 'bursting', 'nonbursting', 0.001)
    set(gcf, 'renderer', 'painters')
    saveas(gcf, fullfile(rootdir, '_analysis', 'R3', 'b_vs_nb.fig'));
    saveas(gcf, fullfile(rootdir, '_analysis', 'R3', 'b_vs_nb.eps'));
    saveas(gcf, fullfile(rootdir, '_analysis', 'R3', 'b_vs_nb.jpg'));
    close(gcf)
    
    [~, r_cellinx, ~] = intersect(vpcells, rhythmic_cells); %FIGR5.B
    [~, nr_cellinx, ~] = intersect(vpcells, nonrhythmic_cells);
    boxstat(BFR(r_cellinx), BFR(nr_cellinx), 'rhythmic', 'non-rhythmic', 0.001)
    set(gcf, 'renderer', 'painters')
    saveas(gcf, fullfile(rootdir, '_analysis', 'R3', 'r_vs_nr.fig'));
    saveas(gcf, fullfile(rootdir, '_analysis', 'R3', 'r_vs_nr.eps'));
    saveas(gcf, fullfile(rootdir, '_analysis', 'R3', 'r_vs_nr.jpg'));
    close(gcf)
    
    [~, s_cellinx, ~] = intersect(vpcells, sync_cellids); %FIGR5.C
    [~, as_cellinx, ~] = intersect(vpcells, nonsync_cellids);
    boxstat(BFR(s_cellinx), BFR(as_cellinx), 'synchronous', 'asynchronous', 0.001)
    set(gcf, 'renderer', 'painters')
    saveas(gcf, fullfile(rootdir, '_analysis', 'R3', 's_vs_as.fig'));
    saveas(gcf, fullfile(rootdir, '_analysis', 'R3', 's_vs_as.eps'));
    saveas(gcf, fullfile(rootdir, '_analysis', 'R3', 's_vs_as.jpg'));
    close(gcf)
end

if figR8spec
    choosecb('VP_tagged_CB')   % choose CellBase with tagged neurons
    resdir_tagged = getpref('cellbase', 'datapath');
    
    response_resdir = fullfile(resdir_tagged,'vpresponsesorter');   % results directory
    acgdir_tagged = fullfile(resdir_tagged,'acg_tagged');
    tagged_vpcells = {'HDB25_180421a_2.1','HDB25_180420a_2.2','HDB26_180524a_4.3','HDB26_180524a_4.2',...
        'HDB26_180524a_4.1','HDB26_180524a_1.1','HDB26_180524a_1.2'}; % cellid list of tagged VP cells
    
    vpacg(tagged_vpcells,acgdir_tagged,true); % autocorrelograms
    
    resdir_tagged_vp = fullfile(rootdir, 'R8', 'response_profiles');
    for k = 1:length(tagged_vpcells)% PSTHs aligned to reward and punishment
        G = figure;
        pause(0.01)
        viewcell2b(tagged_vpcells(k),'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#Feedback','window',[-3 3])
        maximize_figure(G)
        
        cellidt = tagged_vpcells{k};
        cellidt(cellidt=='.') = '_';
        fnm = fullfile(resdir_tagged_vp,[cellidt '_Feedback.jpg']);   % save
        fnm2 = fullfile(resdir_tagged_vp,[cellidt '_Feedback.eps']);   % save
        set(G, 'renderer', 'painters');
        saveas(G,fnm)
        saveas(G,fnm2)
        close(G)
    end
    choosecb('VP_CellBase_B')
end

if figR10spec
    resdirR10 = fullfile(rootdir, '_analysis', 'R10');
    
    if ~isfolder(resdirR10)
        mkdir(resdirR10)
    end
    
    VPsubarea = getvalue('SubArea', vpcells); %subarea information (VPm or VPl)
    
    % FIGR10.B-C - VP e-types in VPvm and VPl
    vpmcells = vpcells(strcmp(VPsubarea, 'VPm')); % proportion of VPvm vs. VPl neurons
    vplcells = vpcells(strcmp(VPsubarea, 'VPl'));
    
    % VPvm e-types
    vpvm_b_cells = intersect(vpmcells, bursting_cells); % bursting VPvm cells
    vpvm_nb_cells = intersect(vpmcells, nonbursting_cells); % non-bursting VPvm cells
    vpvm_r_cells = intersect(vpmcells, rhythmic_cells); % rhythmic VPvm cells
    vpvm_nr_cells = intersect(vpmcells, nonrhythmic_cells); % non-rhythmic VPvm cells
    vpvm_s_cells = intersect(vpmcells, sync_cellids); % synchronous VPvm cells
    vpvm_as_cells = intersect(vpmcells, nonsync_cellids); % asynchronous VPvm cells
    
    % VPl e-types
    vpl_b_cells = intersect(vplcells, bursting_cells); % bursting VPl cells
    vpl_nb_cells = intersect(vplcells, nonbursting_cells); % non-bursting VPl vells
    vpl_r_cells = intersect(vplcells, rhythmic_cells); % rhythmic VPl cells
    vpl_nr_cells = intersect(vplcells, nonrhythmic_cells); % non-rhythmic VPl cells
    vpl_s_cells = intersect(vplcells, sync_cellids); % synchronous VPl cells
    vpl_as_cells = intersect(vplcells, nonsync_cellids); % asynchronous VPl cells
    
    % chi-square statistics
    [~,~,pval_r10_1] = chi2_vp(length(vpvm_b_cells),length([vpvm_b_cells vpvm_nb_cells]),length(vpl_b_cells),length([vpl_b_cells vpl_nb_cells]));
    [~,~,pval_r10_2] = chi2_vp(length(vpvm_r_cells),length([vpvm_r_cells vpvm_nr_cells]),length(vpl_r_cells),length([vpl_r_cells, vpl_nr_cells]));
    [~,~,pval_r10_3] = chi2_vp(length(vpvm_s_cells),length([vpvm_s_cells vpvm_as_cells]),length(vpl_s_cells),length([vpl_s_cells vpl_as_cells]));
    
    
    % FIGR10.D - cue, reward and punishment response of bursting and non-bursting neurons of VPvm and VPl
    vpvm_b_resp = getcellresp(vpvm_b_cells); % bursting VPvm cells
    vpvm_nb_resp = getcellresp(vpvm_nb_cells); % non-bursting VPvm cells
    vpl_b_resp = getcellresp(vpl_b_cells); % bursting VPl cells
    vpl_nb_resp = getcellresp(vpl_nb_cells); % non-bursting VPl cells
    
    % Corresponding chi-square statistics
    [~,~,pval_r10_4] = chi2_vp(length([vpvm_b_resp.cue.excitation' vpvm_b_resp.cue.inhibition']),length(vpvm_b_cells),length([vpvm_nb_resp.cue.excitation' vpvm_nb_resp.cue.inhibition']),length(vpvm_nb_cells));
    [~,~,pval_r10_5] = chi2_vp(length([vpvm_b_resp.reward.excitation' vpvm_b_resp.reward.inhibition']),length(vpvm_b_cells),length([vpvm_nb_resp.reward.excitation' vpvm_nb_resp.reward.inhibition']),length(vpvm_nb_cells));
    [~,~,pval_r10_6] = chi2_vp(length([vpvm_b_resp.punishment.excitation' vpvm_b_resp.punishment.inhibition']),length(vpvm_b_cells),length([vpvm_nb_resp.punishment.excitation' vpvm_nb_resp.punishment.inhibition']),length(vpvm_nb_cells));
    
    [~,~,pval_r10_7] = chi2_vp(length([vpl_b_resp.cue.excitation' vpl_b_resp.cue.inhibition']),length(vpl_b_cells),length([vpl_nb_resp.cue.excitation' vpl_nb_resp.cue.inhibition']),length(vpl_nb_cells));
    [~,~,pval_r10_8] = chi2_vp(length([vpl_b_resp.reward.excitation' vpl_b_resp.reward.inhibition']),length(vpl_b_cells),length([vpl_nb_resp.reward.excitation' vpl_nb_resp.reward.inhibition']),length(vpl_nb_cells));
    [~,~,pval_r10_9] = chi2_vp(length([vpl_b_resp.punishment.excitation' vpl_b_resp.punishment.inhibition']),length(vpl_b_cells),length([vpl_nb_resp.punishment.excitation' vpl_nb_resp.punishment.inhibition']),length(vpl_nb_cells));
    
    % FIGR10.E-G
    % correlation of DV position with rhythmicity indices
    DVpos = getvalue('DVPos', vpcells);
    
    [~, phasic_cellinx, ~] = intersect(vpcells, phasic_cells); % phasic cell indices
    
    regressionplot(BurstIndex, DVpos, 'BurstIndex', 'DVpos', 'regression_DVpos_BI', 'DVpos', resdirR10) % DV position correlation with burst index
    regressionplot(BetaIndex', DVpos(phasic_cellinx), 'BetaIndex', 'DVpos', 'regression_DVpos_BetaI', 'DVpos', resdirR10) % DV position correlation with beta rhythmicity index
    regressionplot(GammaIndex', DVpos(phasic_cellinx), 'GammaIndex', 'DVpos', 'regression_DVpos_GI', 'DVpos', resdirR10) % DV position correlation with gamma rhythmicity index
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% -------------------------------------------------------------------------
function responses = getcellresp(cellids)

% VP neuron responses
cueresp = getvalue('cueresponse',cellids);
rewardresp = getvalue('rewardresponse',cellids);
punishmentresp = getvalue('punishresponse',cellids);
omissionresp = getvalue('omissionresponse',cellids);

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

om_e = cellids(omissionresp == 1);
om_i = cellids(omissionresp == -1);
om_n = cellids(omissionresp == 0);

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
responses.omission.excitation = om_e';
responses.omission.inhibition = om_i';
responses.omission.none = om_n';
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

% -------------------------------------------------------------------------
function regressionplot(xdata, ydata, xtitle, ytitle, filenm, head, resdir)

% Plotting scatter plots
x = xdata;
y = ydata;

X = [ones(length(ydata),1) x];
[b,bint,r,rint,stats] = regress(y,X);
p = stats(3);
pR = corrcoef(y,X(:,2));
R = pR(3);
[b,stats] = robustfit(x,y);

% Regression plot
figure
scatter(x,y);
xlabel(xtitle)
ylabel(ytitle)
title(head)
% setmyplot_Balazs
axis square
icp = b(1);   % intercept
gr = b(2);   % gradient
xx = min(x):0.01:max(x);
yy = xx .* gr + icp;
hold on
plot(xx,yy,'Color',[0.6627 0.6196 0.4039],'LineWidth',2)   % overlay regression line
text('Units','normalized','Position',[0.7 0.7],...
    'String',{['p = ' num2str(p)] ['R = ' num2str(R)]})
set(gcf, 'renderer', 'painters')
% savename = fullfile(resdir, [filenm '.fig']);
savename2 = fullfile(resdir, [filenm '.eps']);
savename3 = fullfile(resdir, [filenm '.jpg']);
% saveas(gcf, savename)
saveas(gcf, savename2)
saveas(gcf, savename3)
end