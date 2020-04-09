function most_pairs = select_monosyn_pairs(filesource)
%	SELECT_MONOSYN_PAIRS finds putative monosynaptic connections within a
%	recording session.
%
%   MOST_PAIRS = SELECT_MONOSYN_PAIRS(FILESOURCE) finds a network of
%   putative monosynaptically connected VP neurons (MOST_PAIRS).
%   Monosynaptic direct connections and synchronous activation was detected
%   based on their cross-correlograms made by VPCCG. Cross correlograms of
%   neuron pairs are copied to RESDIR.

%   See also VPCCG and CCG_GROUPING_VP.

%   Panna Hegedus
%   Institute of Experimental Medicine
%   panna.hegedus@koki.mta.hu
%   05-Feb-2020

%   Code review: BH 2/12/20

% Load_CCG_pairs
load(fullfile(filesource,'cellgroups_tetrodepairs.mat'));
load(fullfile(filesource, 'cellgroups_nontetrodepairs.mat'));

resdir = fullfile(filesource, 'acg_network');
sourcedir = fullfile(filesource, 'allpairs');

if ~isdir(resdir)
    mkdir(resdir)
end

% Find the session with most pairs
Pairs = [monosyn_exc_nttp monosyn_exc_ttp sync_exc_ttp sync_exc_nttp sync_monosyn_nttp]; % pair of cells
numPairs = length(Pairs); % number of Pairs
animalandsession = [];
maxPairsperAnimal = 0;
for i = 1:numPairs
    numsess = 0;
    currentSession = char(Pairs{i}{1}(1:13));
    if sum(strcmp(animalandsession, currentSession)) == 0
        animalandsession = [animalandsession; currentSession ];
        for j=1:numPairs
            currentSession2 = char(Pairs{j}{1}(1:13));
            if sum(strcmp(currentSession, currentSession2)) == 1
                numsess = numsess+1;
            end
        end
        if numsess > maxPairsperAnimal
            maxPairsperAnimal = numsess;
            most_pairs = currentSession;
        end
    end
end

% Copy ccg files into an other folder for easier comparison
for k = 1:numPairs
    currentSession = char(Pairs{k}{1}(1:13));
    if sum(strcmp(most_pairs, currentSession)) == 1
        cellid1 = char(Pairs{k}{1});
        cellid2 = char(Pairs{k}{2});
        tt1=cellid1(15);
        tt2=cellid2(15);
        u1=cellid1(17);
        u2=cellid2(17);
        fnm = ['CCG_' currentSession '_' tt1 '_' u1 '_' currentSession '_' tt2 '_' u2 '.jpg' ];
        copyfile ([sourcedir '\' fnm], [resdir '\' fnm]);
    end
end