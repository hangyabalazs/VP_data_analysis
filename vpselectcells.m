function vpcells = vpselectcells(vpcells)
%VPSELECTCELLS   Select VP neurons from CellBase.
%   VPSELECTCELLS(CELLID) selects cells that are assigned to VP area
%   ('Area1' property set to 'VP') from CELLIDS (default, all cell IDs in
%   CellBase). Only well-isolated units (ID > 20, L_ratio < 0.15) are
%   returned.
%
%   See also SLOWANALYSIS_PAVLOVIAN.

%   Balazs Hangya and Panna Hegedus
%   Laboratory of Systems Neuroscience
%   Institute of Experimental Medicine, Budapest, Hungary

%   Code review: BH 7/24/19, 12/2/19

% Input arguments
narginchk(0,1)
if nargin < 1
    vpcells = [];
end

% Load CellBase
load(getpref('cellbase','fname'));

% List of cellIDs
if isempty(vpcells)
    Lratio = getvalue('Lr_PC');
    ID = getvalue('ID_PC');
    vldty = getvalue('validity');
    area1 = getvalue('Area1');
    isvp = strcmp(area1,{'VP'});   % select VP cells
    ptinx = isvp & ID > 20 & Lratio < 0.15;   % good clusters from VP
    vpcells = CELLIDLIST(ptinx);
else
    if isnumeric(vpcells)  % 'vpcells' can be index set or list of cell IDs
        vpcells = CELLIDLIST(vpcells);
    else
        if ischar(vpcells)
            vpcells = {vpcells};   % only one cell ID
        end
    end
end
vpcells = vpcells(:)';   % convert to row vector