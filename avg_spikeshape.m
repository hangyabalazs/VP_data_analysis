function avg_spikeshape(cellids, resdir, tag)
%AVG_SPIKESHAPE   Plots average waveform of a cell population.
%   AVG_SPIKESHAPE(CELLIDS, RESDIR, TAG) plots average waveform of a cell
%   population listed in CELLIDS. Waveform data is generated by
%   SPIKESHAPE_ANALYSIS.M and added to The Matrix. Result plots are saved
%   to RESDIR with a nametag specifying the cell population (TAG).
%
%   See also SPIKESHAPE_ANALYSIS and ADDANALYSIS

%   Balazs Hangya, Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu
%   18-02-2021

if ~isfolder(resdir) 
    mkdir(resdir)% make results directory
end

spikeshape_data = getvalue('SpikeShape', cellids); % Retrieve data from CellBase

waveform = cell(length(cellids), 1);
waveform = [cellfun(@(x) x{1,1}.Spike, spikeshape_data, 'UniformOutput', false)];
waveform_matrix = cell2mat(waveform);
wave_avg = mean(waveform_matrix);
wave_sd = std(waveform_matrix);

%plot average waveform
H = figure;
errorshade([1:size(waveform_matrix,2)],wave_avg,wave_sd,...
    'LineColor',[0 0 1],'ShadeColor',[0 0 1]) 

set(gcf, 'renderer', 'painters'); % set renderer
ylim([-60 100]); % use uniform Y axis

% save
saveas(gcf, fullfile(resdir, [tag '_avg_spikeshape.eps']));
saveas(gcf, fullfile(resdir, [tag '_avg_spikeshape.jpg']));