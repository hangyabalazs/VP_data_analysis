function [besiinx_cell, biseinx_cell] = redo_vpresponsesorter(cellids, cond, pval, resdir)
%REDO_VPRESPONSESORTER Cell sorting based on burst and single spike response to behavioral events.
%	REDO_VPRESPONSESORTER(CELLIDS, COND, PVAL, RESDIR) is sorting CELLIDS
%	into different response groups based on a given significance level
%	(PVAL) for response magnitude calculated previously by
%	VPRESPONSESORTER_BURST and VPRESPONSESORTER_SINGLE. PSTHs are aligned
%	to behavioral events defined by COND. Output variables BESIINX_CELL and
%	BISEINX_CELL are lists of cellids containing the cells which respond
%	with burst inhibition and single spike activation or vice versa to
%	COND. Cell proportions are plotted on a pie chart and saved to RESDIR.
%
%   See also PSTH_STATS, VPRESPONSESORTER, VPRESPONSESORTER_BURST and
%   VPRESPONSESORTER_SINGLE.

%   Panna Hegedus
%   20-Sept-2020
%   panna.hegedus@koki.mta.hu

%   Code review: BH 10/16/20

% Source directories
source1 = fullfile(getpref('cellbase', 'datapath'),'_paper_figs\code_review3\fig6\PSTH_burst_vs_singlespike\vpresponsesorter_burst');
source2 = fullfile(getpref('cellbase', 'datapath'),'_paper_figs\code_review3\fig6\PSTH_burst_vs_singlespike\vpresponsesorter_single');

switch cond
    case 'cueresponse'
        tag1 = '_StimulusOn_TrialType'; % fing tags for loading psth_stats file
        tag2 = '_StimulusOn_TrialType_single';
    case 'rewardresponse'
        tag1 = '_DeliverAllFeedback_RewardedTrials';
        tag2 = '_DeliverAllFeedback_RewardedTrials_single';
    case 'punishresponse'
        tag1 = '_DeliverAllFeedback_PunishedTrials';
        tag2 = '_DeliverAllFeedback_PunishedTrials_single';
end

[bothiinx, botheinx, besiinx, biseinx, others] = deal(zeros(1)); % preallocate space
besiinx_cell = [];
biseinx_cell = [];
for i = 1:length(cellids)
    cell = cellids{i}; % get cellid
    cell(cell=='.') = '_';
    if exist(fullfile(source1, [cell tag1 '.mat']))
        burst_stat = load(fullfile(source1, [cell tag1 '.mat'])); % load psth_stats
        single_stat = load(fullfile(source2, [cell tag2 '.mat']));
        close all
        
        % Count neuronal responses
        if length(burst_stat.stats1) > 1  % if there are two sounds
            if (burst_stat.stats1{1}.Wpi < pval || burst_stat.stats1{2}.Wpi < pval) && (single_stat.stats1{1}.Wpi < pval || single_stat.stats1{2}.Wpi < pval)
                bothiinx = bothiinx + 1;   % burst and single spike inhibited
            elseif (burst_stat.stats1{1}.Wpa < pval || burst_stat.stats1{2}.Wpa < pval) && (single_stat.stats1{1}.Wpa < pval || single_stat.stats1{2}.Wpa < pval)
                botheinx = botheinx + 1;   % burst and single spike activated
            elseif (burst_stat.stats1{1}.Wpa < pval || burst_stat.stats1{2}.Wpa < pval) && (single_stat.stats1{1}.Wpi < pval || single_stat.stats1{2}.Wpi < pval)
                besiinx = besiinx + 1;   % burst activated single spike inhibited
                besiinx_cell = [besiinx_cell {cellids{i}}];
            elseif (burst_stat.stats1{1}.Wpi < pval || burst_stat.stats1{2}.Wpi < pval) && (single_stat.stats1{1}.Wpa < pval || single_stat.stats1{2}.Wpa < pval)
                biseinx = biseinx + 1;   % burst inhibited and single spike activated
                biseinx_cell = [biseinx_cell {cellids{i}}];
            else
                others = others + 1; % others
            end
        else   % if there was one cue tone
            if (burst_stat.stats1{1}.Wpi < pval ) && (single_stat.stats1{1}.Wpi < pval )
                bothiinx = bothiinx+1;   % burst and single spike inhibited
            elseif (burst_stat.stats1{1}.Wpa < pval ) && (single_stat.stats1{1}.Wpa < pval )
                botheinx = botheinx+1;   % burst and single spike activated
            elseif (burst_stat.stats1{1}.Wpa < pval ) && (single_stat.stats1{1}.Wpi < pval )
                besiinx = besiinx+1;   % burst activated single spike inhibited
                besiinx_cell = [besiinx_cell {cellids{i}}];
            elseif (burst_stat.stats1{1}.Wpi < pval ) && (single_stat.stats1{1}.Wpa < pval )
                biseinx = biseinx+1;   % burst inhibited and single spike activated
                biseinx_cell = [biseinx_cell {cellids{i}}];
            else
                others = others+1;
            end
        end
    end
end

% Create pie chart
X = [bothiinx botheinx besiinx biseinx others];
legend = {};