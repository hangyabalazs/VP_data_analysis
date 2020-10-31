function lineplot_with_errorbars_vp(data, labels, resdir)
%LINEPLOT_WITH_ERRORBARS(DATA,LABELS) Line plot with error bars.
%   LINEPLOT_WITH_ERRORBARS plots DATA as a lineplot (median) with indicated error
%   bars and data label (LABELS). The DATA supposed be in the following format:
%   Each cell in DATA indicates an individual measurement. Rows within each
%   measurement are the data points, while columns are containing data
%   points of different conditions (eg. anticipatory response of two different cues)
%   Results are saved in RESDIR.
%
%   See also SE_OF_MEDIAN

%   Code review: BH 5/11/20

%   Panna Hegedus
%   panna.hegedus@koki.mta.hu
%   30-04-2020

% Plot
NumLines = length(data); % number of datasets to plot
figure;
x = length(data{1}(1,:)); % number of columns in dataset
xvec = [1:x]; % x values to plot
[errors, medians] = deal(nan(1,x)); % errors and median values of each column within the dataset
for i = 1:NumLines
    cDataset = data{i}; % current dataset to plot
    clabel = labels{i}; % label of dataset
    for j = 1:x
        errors(j) = se_of_median(cDataset(:,j)); % error values
        medians(j) = nanmedian(cDataset(:,j)); % median values
    end
    errorbar(xvec, medians, errors) % plot
    hold on
end
legend(labels) % add legend

% Save plot (fig jpg and eps format)
fnm = fullfile(resdir, 'lineplot_with_errorbars.fig');
fnm2 = fullfile(resdir, 'lineplot_with_errorbars.jpg');
fnm3 = fullfile(resdir, 'lineplot_with_errorbars.eps');

saveas(gcf, fnm)
saveas(gcf, fnm2)
set(gcf, 'Renderer', 'painters');
saveas(gcf, fnm3)
close(gcf)