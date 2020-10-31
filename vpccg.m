function vpccg(cellids,resdir,issave)
%VPCCG   Cross-correlation analysis.
%   VPCCG calculates cross-correlations for VP neurons in 50 ms windows.
%   CCG results are saved for non-tetrode pairs and tetrode pairs
%   separately.
%
%   See also CCG.

%   Balazs Hangya, Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu
%   11-Nov-2018

%   Code review: BH 2/12/20, 4/8/20

% Directories
resdir_tetrode = fullfile(resdir,'tetrodepairs');   % results directory for tetrode pairs
resdir_nontetrode = fullfile(resdir,'nontetrodepairs');   % results directory for non-tetrode pairs
if ~isfolder(resdir_tetrode)
    mkdir(resdir_tetrode)
end
if ~isfolder(resdir_nontetrode)
    mkdir(resdir_nontetrode)
end

% Input argument check
narginchk(0,3);
if nargin < 2
    issave = true;   % default saving behavior
end
if nargin < 1
    vpcells = vpselectcells([]);
else
    vpcells = cellids;
end
vpcells = vpcells';

% CCG
segfilter = 'stim_excl_vp';
filterinput = {'light_activation_duration',[-5 5],'margins',[0 0]};
% ccg(vpcells,0.05,'whichcells','nontetrodepairs','resdir',resdir_nontetrode,...
%      'segfilter',segfilter,'filterinput',filterinput,...
%      'minspikeno',100,'maxspikeno',10000,'issave',issave);   % non-tetrode pairs
% ccg(vpcells,0.05,'whichcells','tetrodepairs','resdir',resdir_tetrode,...
%     'segfilter',segfilter,'filterinput',filterinput,...
%     'minspikeno',100,'maxspikeno',10000,'issave',issave);   % tetrode pairs

ccg_grouping_VP(vpcells,'whichcells','nontetrodepairs','resdir',resdir_nontetrode,...
    'segfilter',segfilter,'filterinput',filterinput,...
    'minspikeno',100,'maxspikeno',10000,'issave',issave);   % non-tetrode pairs;
ccg_grouping_VP(vpcells,'whichcells','tetrodepairs','resdir',resdir_tetrode,...
    'segfilter',segfilter,'filterinput',filterinput,...
    'minspikeno',100,'maxspikeno',10000,'issave',issave);   % tetrode pairs;