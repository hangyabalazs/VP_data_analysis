function baseline_vs_burst(vpcells, resdir)
%BASELINE_VS_BURST   Correlates baseline firing rate and bursting.
%   BASELINE_VS_BURST correlates baseline firing rate with burst index and
%   burst spike ratio of neurons listed in VPCELLS. Scatter plots are saved
%   to RESDIR.
%
%   See also EXTRACT_BASELINEFR ACG and VPACG

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   18-02-2021

if ~isfolder(resdir) % make results directory
    mkdir(resdir)
end

NumCells = length(vpcells); % number of cells
baselineFR = getvalue('baseline_FR',vpcells); % baseline firing rate
BI = getvalue('BurstIndex', vpcells); % burst index
propname = 'ISIratio';
if ~ismember(propname,listtag('prop'))
    insertdata([vpcells' num2cell(nan(size(vpcells')))],'type','property','name',propname);
end

ISIratio = nan(1,NumCells);
for i = 1:NumCells
    cellid = vpcells{i};
    
    % Load spike times
    stimes = loadcb(cellid,'Spikes');
    
    % Inter-spike intervals
    isi = diff(stimes);  % inter-spike interval
    isi_100 = isi(isi<0.1);
    isi_100_inx = find(isi<0.1); % reference to ISI<100ms
    stimes_100 = stimes(unique([isi_100_inx isi_100_inx+1]));
    burstinx = find(isi_100<0.01);   % ISI < 10 ms
    
    % Detect bursts
    % Burst: first ISI < 10 ms, subsequent ISIs < 15 ms
    bursts = {};
    used = [];
    burst1 = {};
    for k = 1:length(burstinx)
        if ismember(burstinx(k),used)
            continue   % already included in the previous burst
        end
        bursts{end+1} = stimes_100([burstinx(k), burstinx(k)+1]);
        next = burstinx(k) + 1;
        if next > length(isi_100)
            burst1{end+1} = bursts{end}(1);
            continue   % last ISI of the cell
        end
        while isi_100(next) < 0.015  % if conseq. ISIs < 15 ms
            bursts{end} = [bursts{end}; stimes_100(next+1)];
            used = [used next]; %#ok<AGROW>
            next = next + 1;
            if next > length(isi_100)
                break   % last ISI of the cell, quit while loop
            end
        end
        burst1{end+1} = bursts{end}(1);
    end
    bursts = vertcat(bursts{:});
    single = setdiff(stimes_100, bursts);
    
    % Create ISI ratio
    ISIratio(i) = length(bursts) / length(stimes_100); % burst spikes / all spikes
    st = setvalue(cellid,propname,ISIratio(i));   % insert ISI ratio into CB
end

% Make scatter plots
figure; % Baseline firing rate correlated with burst spike ratio
scatter(baselineFR(BI>0.2), ISIratio(BI>0.2), 'bo'); % bursting cells
hold on
scatter(baselineFR(BI<=0.2), ISIratio(BI<=0.2), 'ro'); % nonbursting cells
set(gcf, 'renderer', 'painters');
saveas(gcf, fullfile(resdir, 'BFR_vs_ISI.fig'));
saveas(gcf, fullfile(resdir, 'BFR_vs_ISI.eps'));
saveas(gcf, fullfile(resdir, 'BFR_vs_ISI.jpg'));
close(gcf)

figure; % Baseline firing rate correlated with burst index
scatter(baselineFR(BI>0.2), BI(BI>0.2), 'bo'); % bursting cells
hold on
scatter(baselineFR(BI<=0.2), BI(BI<=0.2), 'ro'); % nonbursting cells

set(gcf, 'renderer', 'painters');
saveas(gcf, fullfile(resdir, 'BFR_vs_BI.fig'));
saveas(gcf, fullfile(resdir, 'BFR_vs_BI.eps'));
saveas(gcf, fullfile(resdir, 'BFR_vs_BI.jpg'));
close(gcf)