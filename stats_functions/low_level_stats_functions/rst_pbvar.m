function [diff,CI,p,sig] = rst_pbvar(data1,data2)
%[diff,CI,P,SIG] = pbvar(data1,data2)
%PBVAR - percentile bootstrap test to compare the variances of 2
%independent groups
% 2 tailed bootstrap test based on the estimation of H1, the
% hypothesis of an experimental effect. nboot the number of bootstrap
% samples = 600. 
% The endpoints of the CI are computed using Wilcox's
% adjustements for small sample sizes to provide a better control over the
% probability of a type I error. 
% Alpha = .05
% mu = 0 - null hypothesis
%
% INPUTS:
% -------
%
% DATA1 and DATA2 are 2 vectors
%
% OUTPUTS:
% --------
%
% LOCI/HICI: low and high boundaries of the confidence interval for
% the difference between 2 groups. Contrary to permtest and
% pbH0, the confidence interval calculated here is not under the
% null hypothesis. Therefore, instead of being centred on zero, it will
% follow the difference between the 2 experimental conditions compared.
% As a simple decisional rule, when the confidence interval under H1 does
% not include zero, the difference is considered significant.
%
% P value of the effect
%
% SIG: significativity of the experimental difference, binary output,
% 1=YES, 0=NO
%
% Guillaume A. Rousselet - University of Glasgow - Dec 2008
%
% See also PBH1

rand('state',sum(100*clock));

% if nargin<3 || isempty(nboot)
%   nboot=1000;
% end
% if nargin<4 || isempty(alpha)
%     alpha=.05;
% end

nboot = 600;
n1 = length(data1);
n2 = length(data2);
% lo = round(nboot.*alpha./2); % get CI boundary indices
% hi = nboot - lo;
nmin = min(n1,n2);
if nmin<40
    lo = 6; hi = 593;
elseif nmin < 80
    lo = 7; hi = 592;
elseif nmin < 180
    lo = 10; hi = 589;
elseif nmin < 250
    lo = 13; hi = 586;
elseif nmin >= 250
    lo = 15; hi = 584;
end
bootdiff = zeros(1,nboot);

for kk = 1:nboot % bootstrap with replacement loop
    bootdiff(kk)=var(randsample(data1,nmin,true)) - var(randsample(data2,nmin,true));
end

diffsort = sort(bootdiff); % sort in ascending order
CI(1) = diffsort(lo+1);
CI(2) = diffsort(hi);

diff = var(data1) - var(data2);

if nargout > 2
    mu = 0;
    alpha = .05;
    p = length(find(bootdiff>mu))./nboot;
    p = 2.*min(p,1-p);
    sig = p<=alpha; 
end

return
            
      

       