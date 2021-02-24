function compare_spikeshape(cellids, labels, resdir)
%COMPARE_SPIKESHAPE   Performing statistical comparison of spike shape
%properties across neuronal populations.
%   COMPARE_SPIKESHAPE(CELLIDS, LABELS, RESDIR) compares spike shape
%   properties calculated and added to CellBase by SPIKESHAPE_ANALYSIS_P.M.
%   CELLIDS is a 1x2 cell containing two lists of cellids. LABELS is a list
%   of string containing the label for the cellid groups defined in
%   CELLIDS. Mann-Whitney U tests are performed comparing each feature. The
%   results are visualized az boxplots and are saved to RESDIR.

%   See also SPIKESHAPE_ANALYSIS_P, SPIKESHAPEANALYSIS and BOXSTAT.

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   18-02-2021

if ~isfolder(resdir) % If results directory not exists, make one
    mkdir(resdir)
end

cell_subset_type = cellids{1}; % two cell groups to compare
cell_subset_nontype = cellids{2};

spikeshape_data_type = getvalue('SpikeShape', cell_subset_type); % Retrieve data from CellBase
spikeshape_data_nontype = getvalue('SpikeShape', cell_subset_nontype);

FieldNames = fieldnames(spikeshape_data_type{1, 1}{1, 1}); % Fieldnames
NumFeatures = length(FieldNames); % number of Features to compare

for i = 2:NumFeatures %loop through features
    cFeature_type = nan(length(spikeshape_data_type),1);
    cFeature_nontype = nan(length(spikeshape_data_nontype),1);
    
    cFeature_name = FieldNames{i}; % choose a property
    
    cFeature_type = [cellfun(@(x) x{1,1}.(cFeature_name), spikeshape_data_type)]; % get the corresponding data for each group
    cFeature_nontype = [cellfun(@(x) x{1,1}.(cFeature_name), spikeshape_data_nontype)];
    
    if isnumeric(cFeature_type) % only compare numeric features (e.g do not perform statistics on waveform data)
    boxstat(cFeature_type, cFeature_nontype, labels{1}, labels{2}, 0.005);     %boxplot and statistics
    title(cFeature_name)
    
    % save
    saveas(gcf, fullfile(resdir, ['comapre_spikefeature_' cFeature_name '.fig']))
    saveas(gcf, fullfile(resdir, ['comapre_spikefeature_' cFeature_name '.jpg']))
    set(gcf, 'renderer', 'painters');
    saveas(gcf, fullfile(resdir, ['comapre_spikefeature_' cFeature_name '.eps']))
    close(gcf)
    end
end