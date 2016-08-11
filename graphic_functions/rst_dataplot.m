function rst_dataplot(Data,options)

% generic function to plot data per groups
% 
% FORMAT rst_dataplot(Data,options)
%
% INPUTS Data is a matrix, and data are taken column wise
%        options are 'scatter', 'summary' and 'distrib'
%        'scatter' is followed by the same arguments as the scatter function
%        'summary' can be 'mean', 'trimmean', 'median'
%                  any summary value comes with 95% High Density Intervals
%                  (derived from a Bayesian bootsrap)
%        'distrib' plot the estimated kernel density using a random shifted
%        histogram approach (i.e. fully non parametric and bounded)

