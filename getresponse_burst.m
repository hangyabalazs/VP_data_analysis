function getresponse_burst(cellids, cond, resdir)
% GETRESPONSE_BURST shows single spike and burst response to cue reward and
% punishment respectively.
% GETRESPONSE_BURST(CELLIDS, COND) finds cells (from CELLIDS) with burst and/or single
% spike response to cue reward and punishment (COND) respectively. Results are
% visualized on a pie chart and saved to RESDIR.

% see also GETVALUE


switch cond
    case 'cue'
        burst = getvalue('cueresponse_burst', cellids);
        single = getvalue('cueresponse_single', cellids);
    case 'reward'
        burst = getvalue('rewardresponse_burst', cellids);
        single = getvalue('rewardresponse_single', cellids);
    case 'punishment'
        burst = getvalue('punishresponse_burst', cellids);
        single = getvalue('punishresponse_single', cellids);
end

botheinx = find(burst ==1 & single ==1);
bothiinx = find(burst == -1 & single == -1);
beinx = find(burst == 1 & single == 0);
biinx = find(burst == -1 & single == 0);
seinx = find(burst == 0 & single == 1);
siinx = find(burst == 0 & single == -1);
besinx = find(burst == 1 & single == -1);
biseinx = find(burst == -1 & single == 1);
n = find(burst == 0 & single == 0);

both_active = length((botheinx)); % both bursts and single spikes are activated
both_inhibit = length((bothiinx)) ; % both bursts and single spikes are inhibited
burst_active = length((beinx));
burst_inh = length((biinx));
single_active = length((seinx));
single_inh = length((siinx));
burst_active_single_inh = length((besinx));
single_active_burst_inh = length((biseinx)); 
none = length(n);

X = [both_active both_inhibit burst_active burst_inh single_active single_inh burst_active_single_inh single_active_burst_inh none];
pie(X)