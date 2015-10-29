function [BEST,EST,CI,p] = onesampb(x,nboot,alpha,est,q,mu)

% function [BEST,EST,CI,p] = onesampb(x,nboot,alpha,est,q,mu)
%
% Compute a percentile bootstrap confidence interval
%
%   INPUTS:
%           x (vector)
%           nboot (resamples, default 1000)
%           alpha (default 5%)
%           est (estimator, default 'median')
%           q (argument for estimator, default .5 for Harrell-Davis
%           estimator, 0.2 for trimmed mean
%           mu (optional null hypothesis, necessary to get a p-value in
%           one-sample hypothesis testing procedure, default 0)
%
%   OUTPUTS:
%           BEST - bootstrap estimates
%           EST - estimate
%           CI - confidence interval
%           p - bootstrap p value
%
% GAR - University of Glasgow - Dec 2007
% added EST output - Nov 2008
% added BEST output - Feb 2012

if nargin<2 || isempty(nboot)
    nboot=1000;
end
if nargin<3 || isempty(alpha)
    alpha=.05;
end
if nargin<4 || isempty(est)
    est='median';
end

if nargin<5 || isempty(q)
    switch lower(est)
        case {'hd'} % Harrell-Davis estimator
            q=.5;
        case {'tm'} % trimmed mean
            q=20;
        case {'trimmean'} % trimmed mean - Matlab function
            q=20; % later on multiplied by 2 to achieve 20% trimming of both tails
    end
end

switch lower(est)
    case {'hd'} % Harrell-Davis estimator
        eval(['EST=',est,'(x,q);'])
    case {'tm'} % trimmed mean
        eval(['EST=',est,'(x,q);'])
    case {'trimmean'} % trimmed mean - Matlab function
        eval(['EST=',est,'(x,q.*2);'])
    otherwise
        eval(['EST=',est,'(x);'])
end

n = length(x);
lo = round(nboot.*alpha./2);
hi = nboot - lo;

for kk = 1:nboot % bootstrap with replacement loop
    switch lower(est)
        case {'hd'} % Harrell-Davis estimator
            eval(['bootx(kk)=',est,'(x(randi(n,n,1)),q);'])
        case {'tm'} % trimmed mean
            eval(['bootx(kk)=',est,'(x(randi(n,n,1)),q);'])
        case {'trimmean'} % trimmed mean - Matlab function
            eval(['bootx(kk)=',est,'(x(randi(n,n,1)),q.*2);'])
        otherwise
            eval(['bootx(kk)=',est,'(x(randi(n,n,1)));'])
    end
end

diffsort = sort(bootx); % sort in ascending order
CI(1) = diffsort(lo+1);
CI(2) = diffsort(hi);
BEST = bootx;

if nargout > 2
    if nargin < 6 || isempty(mu)
        mu=0;
    end
    p = sum(bootx>mu)./nboot + sum(bootx==mu)./(2.*nboot);
    p = 2.*min(p,1-p);
end
