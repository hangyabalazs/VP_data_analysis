function [sync_cellids, nonsync_cellids] = sync_pair_responses(vpcells,sourcedir)
%SYNC_PAIR_RESPONSES   Find synchronous neurons and lists their cellids.
%   [SYNC, NONCSYNC, SYNC_CELLIDS, NONSYNC_CELLIDS] =
%   SYNC_PAIR_RESPONSES(VPCLELS,SOURCEDIR) is sorting cells from VPCELLS
%   into synchronous and nonsyncronous groups (SYNC_CELLIDS,
%   NONSYNC_CELLIDS) based on their cross-correlograms generated by VPCCG.
%
%   See also VPCCG.

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   05-Feb-2020 

%   Code review: BH 4/8/20

% Load CellBase
loadcb
if isempty(vpcells)
    vpcells = vpselectcells([]);
end

% Load_CCG_pairs
load(fullfile(sourcedir,'\cellgroups_tetrodepairs.mat'));
load(fullfile(sourcedir,'\cellgroups_nontetrodepairs.mat'));

% Generate a list of neurons that participate in the pairs
Pairs = [sync_exc_nttp sync_exc_ttp sync_monosyn_nttp];
sync_cellids = unique([Pairs{:}]);

sync_cellids = vpcells(ismember(vpcells,sync_cellids)); % nonsynchronous cells
sync_cellinx = find(ismember(CELLIDLIST ,sync_cellids)); % nonsynchronous cellids

nonsync_cellids = vpcells(~ismember(vpcells,sync_cellids)); % nonsynchronous cells
nonsync_cellinx = find(ismember(CELLIDLIST ,nonsync_cellids)); % nonsynchronous cellids