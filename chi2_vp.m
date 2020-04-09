function [tbl,chi2stat,pval] = chi2_vp(proportion1, group1,proportion2, group2)
%CHI2_VP   Chi square test for VP data.
%	CHI2_VP(GROUP1, GROUP2) is performing chi squre test for homogeneity of
%	GROUP1 and GROUP2 (2x2 contingency table).
%
%   See also CROSSTAB.

% Observed data
n1 = proportion1; 
N1 = group1;
n2 = proportion2; 
N2 = group2;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)
