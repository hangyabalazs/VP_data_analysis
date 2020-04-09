function vpacg(cellids,resdir,issave)
%VPACG   Auto-correlation analysis.
%   VPACG calculates auto-correlations for VP neurons in 500 ms windows.
%   ACG results are saved.
%
%   See also ACG.

%   Balazs Hangya, Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu
%   11-Nov-2018

%   Code review: BH 2/12/20, 4/8/20

% Directories
if ~isfolder(resdir)
    mkdir(resdir)
end

% Input argument check
narginchk(0,2);
if nargin < 2
    issave = true;   % default saving behavior
end
if nargin < 1
    vpcells = vpselectcells([]);
else
    vpcells = cellids;
end

% ACG
segfilter = 'stim_excl_vp';
filterinput = {'light_activation_duration',[-5 5],'margins',[0 0]};
acg(vpcells,0.5,'resdir',resdir,...
     'segfilter',segfilter,'filterinput',filterinput,...
     'minspikeno',100,'maxspikeno',10000,'issave',issave);
