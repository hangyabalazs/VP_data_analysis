function [M, H] = response_matrix(cellids)
%RESPONSE_MATRIX   Number of neurons with all combinations of responses.
%   RESPONSE_MATRIX(CELLIDS) returns the number of cells from CELLIDS with
%   all combinations of positive, negative, or no responses to cue, reward
%   and punishment. Individual responses should be stored in CellBase using
%   1, -1 and 0 codes.
%
%   See also ULTIMATE_PSTH and PSTH_STATS.

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   09-Dec-2019

%   Code review: BH 12/9/19

% Individual responses
cueresp = getvalue('cueresponse',cellids);  % cue
rewardresp = getvalue('rewardresponse',cellids);  % reward
punishmentresp = getvalue('punishresponse',cellids);  % punishment

% Aggregate matrix
M = nan(9,3);
for minx = 1:27
    c = mod(ceil(minx/9)+1,3)-1;  % cue
    r = mod(ceil(minx/3)+1,3)-1;  % reward
    p = mod(minx+1,3)-1;  % punishment
    M(minx) = sum(rewardresp==r&cueresp==c&punishmentresp==p);
end
H = imagesc(M);