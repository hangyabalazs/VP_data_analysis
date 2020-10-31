function sorted_struct = sortresponses_animal(animals, cell_struct)
%SORTRESPONSES_ANIMAL   Get cell responses for cells grouped by animals.
%   SORTED_STRUCT = SORTRESPONSES_ANIMAL(ANIMALS, CELL_STRUCT) gets cell
%   responses for cells grouped by animals (ANIMALS) listed in CELL_STRUCT.
%   Result struct (SORTED_STRUCT) is returned as an output.
%
%   See also SLOWANALYSIS_PAVLOVIAN.

%   Panna Hegedus
%   hegedus.panna@koki.hu
%   30-09-2020

%   Code review: BH 10/16/20

cell_struct_fields = struct2cell(cell_struct);
sorted_struct = struct;
for i = 1:length(cell_struct_fields) % loop through animals
    if ~isempty(cell_struct_fields{i})
        sorted_struct.(animals{i}) = getcellresp(cell_struct_fields{i});
    end
end

% -------------------------------------------------------------------------
function responses = getcellresp(cellids)

% VP neuron responses
cueresp = getvalue('cueresponse',cellids);
rewardresp = getvalue('rewardresponse',cellids);
punishmentresp = getvalue('punishresponse',cellids);

% Cue response
cue_e = cellids(cueresp == 1);  % activated
cue_i = cellids(cueresp == -1); % inhibited
cue_n = cellids(cueresp == 0);  % non-responsive

%reward response
rew_e = cellids(rewardresp == 1);  % activated
rew_i = cellids(rewardresp == -1); % inhibited
rew_n = cellids(rewardresp == 0);  % non-responsive

%punishment response
pun_e = cellids(punishmentresp == 1);  % activated
pun_i = cellids(punishmentresp == -1); % inhibited
pun_n = cellids(punishmentresp == 0);  % non-responsive

% Output
responses = struct;
responses.cue.excitation = cue_e';
responses.cue.inhibition = cue_i';
responses.cue.none = cue_n';
responses.reward.excitation = rew_e';
responses.reward.inhibition = rew_i';
responses.reward.none = rew_n';
responses.punishment.excitation = pun_e';
responses.punishment.inhibition = pun_i';
responses.punishment.none = pun_n';