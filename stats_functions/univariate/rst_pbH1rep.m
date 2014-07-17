function [CI,p,sig,padj] = rst_pbH1rep(data1,data2,nboot,alpha,est,q,mu)

%PBH1REP - percentile bootstrap test to compare 2 dependent groups
% [CI,p,sig,padj] = pbH1rep(data1,data2,nboot,alpha,est,q,mu)
% 2 tailed bootstrap test based on the estimation of H1, the
% hypothesis of an experimental effect.
% The 2 groups must have the same size.
%
% INPUTS:
% -------
%
% DATA1 and DATA2 are 2 vectors:
% 
%           nboot (resamples, default 1000)
%           alpha (default 5%)
%           est (estimator, default 'median')
%           q (argument for estimator, default .5 for Harrell-Davis
%           estimator, 20% for trimmed mean
%           mu (optional null hypothesis, default 0)
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
% PADJ: adjusted p-value
% The percentile bootstrap method above can be unsatisfactory when making
% inferences about means, variances, and in the case of small sample sizes
% with M-estimators or the Harrell-Davis estimate of the median.
% The risk is to have an actual type I error rate lower than the nominal
% level. We can correct for that by computing a bias-adjusted p-value.
% See Wilcox p.192-193
%
% Guillaume A. Rousselet - University of Glasgow - January 2008
% added bias-adjusted p-value - Nov 2008
%
% See also PBH1 HD TM TRIMMEAN

rand('state',sum(100*clock));

if nargin<3 || isempty(nboot)
  nboot=1000;
end
if nargin<4 || isempty(alpha)
    alpha=.05;
end
if nargin<5 || isempty(est)
est='median';
end

if nargin<6 || isempty(q)
    switch lower(est)
        case {'hd'} % Harrell-Davis estimator
            q=.5;
        case {'tm'} % trimmed mean
            q=20;
        case {'trimmean'} % trimmed mean - Matlab function
            q=20; % later on multiplied by 2 to achieve 20% trimming of both tails
        otherwise
            q=0;
    end
end

n1 = length(data1);
n2 = length(data2);
lo = round(nboot.*alpha./2); % get CI boundary indices
hi = nboot - lo;

for kk = 1:nboot % bootstrap with replacement loop
    bootsamp = randsample(1:n1,n1,true);
    switch lower(est)
    case {'hd'} % Harrell-Davis estimator
            eval(['bootdiff(kk)=',est,'(data1(bootsamp),q) -',est,'(data2(bootsamp),q);'])
    case {'tm'} % trimmed mean
            eval(['bootdiff(kk)=',est,'(data1(bootsamp),q) -',est,'(data2(bootsamp),q);'])
    case {'trimmean'} % trimmed mean - Matlab function
            eval(['bootdiff(kk)=',est,'(data1(bootsamp),q.*2) -',est,'(data2(bootsamp),q.*2);'])
    otherwise
            eval(['bootdiff(kk)=',est,'(data1(bootsamp)) -',est,'(data2(bootsamp));'])
    end    
end

diffsort = sort(bootdiff); % sort in ascending order
CI(1) = diffsort(lo+1);
CI(2) = diffsort(hi);

if nargout > 1 % Compute P values
    if nargin < 7 || isempty(mu)
        mu=0;
    end
    p = length(find(bootdiff>mu))./nboot;
    p = 2.*min(p,1-p);
    sig = p<=alpha;

switch lower(est)
    case {'hd'} % Harrell-Davis estimator
        eval(['cdata1=data1-',est,'(data1,q);'])
        eval(['cdata2=data2-',est,'(data2,q);'])
    case {'tm'} % trimmed mean
        eval(['cdata1=data1-',est,'(data1,q);'])
        eval(['cdata2=data2-',est,'(data2,q);'])
    case {'trimmean'} % trimmed mean - Matlab function
        eval(['cdata1=data1-',est,'(data1,q.*2);'])
        eval(['cdata2=data2-',est,'(data2,q.*2);'])
    otherwise
        eval(['cdata1=data1-',est,'(data1);'])
        eval(['cdata2=data2-',est,'(data2);'])
end

for kk = 1:nboot % bootstrap with replacement loop
    bootsamp = randsample(1:n1,n1,true);
    switch lower(est)
        case {'hd'} % Harrell-Davis estimator
            eval(['bootdiff(kk)=',est,'(cdata1(bootsamp),q) -',est,'(cdata2(bootsamp),q);'])
        case {'tm'} % trimmed mean
            eval(['bootdiff(kk)=',est,'(cdata1(bootsamp),q) -',est,'(cdata2(bootsamp),q);'])
        case {'trimmean'} % trimmed mean - Matlab function
            eval(['bootdiff(kk)=',est,'(cdata1(bootsamp),q.*2) -',est,'(cdata2(bootsamp),q.*2);'])
        otherwise
            eval(['bootdiff(kk)=',est,'(cdata1(bootsamp)) -',est,'(cdata2(bootsamp));'])
    end
end

qstar = length(find(bootdiff>0))./nboot;
qstar = 2.*min(qstar,1-qstar);

% bias-adjusted p-value:
padj = p-.1.*(qstar-.5); % with increasing sample size the adjustment becomes negligeable
padj = 2.*min(padj,1-padj);

end

return




            
      

       