function vp_spikeshape_analysis(resdir)
%VP_SPIKESHAPE_ANALYSIS   Performs spike shape analysis on all cells of CELLIDLIST.
%   VP_SPIKESHAPE_ANALYSIS(RESDIR) performs spike shape analysis on all cells of
%   CELLIDLIST. Results are added to The Matrix and are saved to RESDIR.
%
%   See also SPIKESHAPE_ANALYSIS, SPIKESHAPE_ANALYSIS_P and ADDANALYSIS.

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   18-02-2021

if ~isfolder(resdir) % make results directory
    mkdir(resdir)
end

% Perform analysis and add waveform features to CellBase
addanalysis(@spikeshapeanalysis_p,...
  'mandatory',{true, resdir},...
  'property_names',{'SpikeShape_test'});
end