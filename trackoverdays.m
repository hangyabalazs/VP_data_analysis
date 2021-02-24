function Cells2Exclude = trackoverdays(cellids)
%TRACKOVERDAYS   Determine if the same neuron was recorded in different sessions.
%   TRACKOVERDAYS(CELLIDS) estimates the probability that the two clusters,
%   recorded on different days but on the same tetrode, with 150 microns 
%   based on track reconstruction, belong to the same neuron. The method is
%   based on Fraser and Schwartz (2012). Three similarity scores are used:
%   maximum waveform crosscorrelation (Fisher-transformed), autocorrelogram
%   correlation (Fisher-transformed) and absolute log difference of firing
%   rates. The scores are tested against a bootstrap distribution of pairs
%   of cells recorded in different animals. p=0.05 is used for the first
%   two and p=0.1 is used for the third scores (pairs have to pass all
%   three to be considered the same neuron). The firing rate criterion is
%   more laxed than the other two, because firing rates can change between
%   and even within sessions. The crosscorrelation criterion was omitted
%   from the Fraser paper, since neighbors can change when the electrodes
%   are moved.
%
%   Reference:
%   Fraser GW, Schwartz AB (2012) Recording from the same neurons 
%   chronically in motor cortex, J Neurophys, 107: 1970 –1978.

%   Balazs Hangya, 4-Feb-2021
%   Institute of Experimental Medicine
%   hangya.balazs@koki.hu 

% Input argument check
if nargin < 1
    cellids = vpselectcells([]);
end

% Initialize bootstrap null distribution
NumCells = numel(cellids);  % number of neurons
prms = nchoosek(1:NumCells,2);  % all possible pairs
rprms = prms(randperm(size(prms,1)),:);  % randomize order
mx = 10000;   % maximum of pairs considered, to save CPU time
c1 = cellids(rprms(1:mx,1));  % candidates for the first neurons of pairs
c2 = cellids(rprms(1:mx,2));  % candidates for the second neurons of pairs
r1 = getvalue('RatId',c1);  % animal ID of the first neurons
r2 = getvalue('RatId',c2);  % animal ID of the second neurons
diffmouseinx = find(r1-r2~=0);  % when the pair is NOT from the same mouse

bno = 1000;  % bootstrap sample size
PairsOfCells = [c1(diffmouseinx(1:bno))'...
    c2(diffmouseinx(1:bno))'];   % take a bootstrap sample; pairs are NOT from the same mouse

Bootstrap null distribution
[FWR, FAR, FFR] = deal(nan(1,bno));
for iP = 1:bno
    disp(iP/bno)
    cellid1 = PairsOfCells{iP,1};  % first cellid in the pair
    cellid2 = PairsOfCells{iP,2};  % second cellid in the pair
    [FWR(iP), FAR(iP), FFR(iP)] = simscores(cellid1,cellid2);   % Fisher-transformed waveform and ACG correlation
    close all
end
% keyboard
% load('simscore_distributions3.mat')

% Critical values
cFWR = prctile(FWR,95);   % critical value corresponding to the upper 0.05
cFAR = prctile(FAR,95);
cFFR = prctile(FFR,10);

% Potential duplicates
[Duplicates, Tested_Pairs] = deal({});
Tested_Pairs_Scores = [];
for iC = 1:NumCells
    disp(iC)
    cellid = cellids{iC};
    % [animalID, sessionID, Tetrode, Unit] = cellid2tags(cellid);
    RatID = getvalue('RatId',cellid);  % animal ID
    datenum = getvalue('DateNum',cellid);  % recording date
    tetrode = getvalue('Tetrode',cellid);  % recording tetrode
    DV = getvalue('DVpos',cellid);  % recording depth
    alldatenum = getvalue('DateNum',cellids);  % recording date for all cells
    allRatIDs = getvalue('RatId',cellids);  % animal ID for all cells
    alltetrode = getvalue('Tetrode',cellids);  % recording tetrode for all cells
    allDV = getvalue('DVpos',cellids);
    samemouseinx = allRatIDs == RatID;  % indices of cells from the same mouse
    diffsessioninx = alldatenum > datenum;  % indices of cells from later sessions
    sametetrodeinx = tetrode == alltetrode;  % indices of cells recorded on the same tetrode
    potential_duplicates = cellids(samemouseinx&sametetrodeinx&...
        diffsessioninx&allDV-DV>=0&allDV-DV<=150);  % cell in the same mouse, on the same tetrode but different session,
    % within 150 um but not above the neuron
    
    % Test the similarity scores of potential duplicates
    NumDup = length(potential_duplicates);  % number of potential duplicates
    for iD = 1:NumDup   % loop through potential duplicates
        cellidd = potential_duplicates{iD};
        [dFWR, dFAR, dFFR] = simscores(cellid,cellidd);   % Fisher-transformed waveform and ACG correlation
        close all
        Tested_Pairs = [Tested_Pairs; {cellid cellidd}]; %#ok<AGROW>  % store similarity score values
        Tested_Pairs_Scores = [Tested_Pairs_Scores; dFWR, dFAR, dFFR]; %#ok<AGROW>  % store similarity score values
        if dFWR > cFWR && dFAR > cFAR && dFFR < cFFR  % if all similarity scores are above critical value (below for FR similarity)
            Duplicates = [Duplicates; {cellid cellidd}];   %#ok<AGROW> % add to list of 'duplicates'
        end
        save('tested_pairs','Tested_Pairs','Tested_Pairs_Scores')
    end
end

keyboard

% delinx = zeros(1,size(Duplicates,1));
% for k = 1:size(Duplicates,1)
%     datenum1 = getvalue('DateNum',Duplicates(k,1));  % recording date for the first cell
%     datenum2 = getvalue('DateNum',Duplicates(k,2));  % recording date for the second cell
%     if datenum2 < datenum1  % pairs may be counted twice
%         delinx(k) = 1;
%     end
% end
% Duplicates(logical(delinx),:) = [];   % remove duplicated duplicates

% delinx = zeros(1,size(Duplicates,1));
% for k = 1:size(Duplicates,1)
%     fr1 = getvalue('baseline_FR',Duplicates(k,1));  % FR for the first cell
%     fr2 = getvalue('baseline_FR',Duplicates(k,2));  % FR for the second cell
%     dFFR = abs(log(fr1)-log(fr2));  % log difference (assumed log-normal FR distribution)
%     if dFFR > cFFR  % drop pairs with FR difference above 10% bootstrap critical value
%         delinx(k) = 1;
%     end
% end
% Duplicates(logical(delinx),:) = [];   % remove cellx with large FR difference

Cells2Exclude = unique(Duplicates(:,2));

% -------------------------------------------------------------------------
function [FWR, FAR, FFR] = simscores(cellid1,cellid2)

% Waveform correlation
wave1 = extractSpikeWaveforms(cellid1,'all','chans','mean_all');   % average waveform on all channels
wave2 = extractSpikeWaveforms(cellid2,'all','chans','mean_all');
twave1 = wave1';  % transpose
twave2 = wave2';
catwave1 = twave1(:)';  % concatenate channels
catwave2 = twave2(:)';
ncatwave1 = catwave1 ./ max(catwave1);   % normalize to a maximum of 1 (remove the effect of proportional amplitude changes)
ncatwave2 = catwave2 ./ max(catwave2);
[xr, lags] = xcorr(ncatwave1,ncatwave2);  % cross-correlogram of waveforms
hv = floor(size(wave1,2)/2);  % half of the waveform 'length' (about half ms)
WR = max(xr(lags>=-hv&lags<=hv));  % waveform correlation: maximal cross-correlation within half-wave times from 0 lag
WR = WR / length(xr);  % normalize between -1 and 1
FWR = atanh(WR);  % Fisher's Z-transform (approximately normal distribution)

% Autocorrelogram
ac1 = acg(cellid1,0.1,'dt',0.005,'minspikeno',0,'maxspikeno',100000);  % autocorrelogram: 100 ms window, 5 ms resolution
[ac2, aclags] = acg(cellid2,0.1,'dt',0.005,'minspikeno',0,'maxspikeno',100000);
acc = corrcoef(ac1(aclags>0),ac2(aclags>0));  % Pearson's correlation
AR = acc(1,2);  % ACG correlation
FAR = atanh(AR);  % Fisher's Z-transform (approximately normal distribution)

% Baseline firing rate
fr1 = getvalue('baseline_FR',cellid1);  % baseline firing rate
fr2 = getvalue('baseline_FR',cellid2);
FFR = abs(log(fr1)-log(fr2));  % log difference (assumed log-normal FR distribution)