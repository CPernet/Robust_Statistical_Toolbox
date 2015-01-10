function [pID,pN] = rst_FDR(p,alpha)
% 
% p   - vector of p-values
% alpha   - False Discovery Rate level
%
% pID - p-value threshold based on independence or positive dependence
% pN  - "Nonparametric" (any covariance structure) p-value threshold
%______________________________________________________________________________
% Based on FDR.m     1.4 Tom Nichols 02/07/02
% short coded by C. Pernet.

p(isnan(p)) = [];
p = sort(p(:));
V = length(p);
I = (1:V)';

cVID = 1;
cVN = sum(1./(1:V));

pID = p(max(find(p<=I/V*alpha/cVID)));
pN  = p(max(find(p<=I/V*alpha/cVN)));

return