function [EST,CI,bootx,p] = rst_pbCI(X,nboot,alphav,est,q,mu)

% function [EST,CI,bootx,p] = rst_pbCI(x,nboot,alpha,est,q,mu)
% 
% Compute a percentile bootstrap confidence interval
%   INPUTS:
%           X is a matrix of data
%           nboot number of resamples (default 1000)
%           alphav is the alpha level (default 5% bonferroni adjustred for the number of columns of x)
%           est is the choosen estimator: 
%               - hd (Harrell-Davis, default q=0.5 ie HD of the median) 
%               - tm (default 20%, same as trimmean but using our function) 
%               - trimmen (default 20%, same as tm but using matlab function)
%               - median (default)
%           q argument for estimator
%           mu (optional null hypothesis, necessary to get a p-value in
%           one-sample hypothesis testing procedure, default 0)
%
%   OUTPUTs:
%          EST is the estimator value
%          CI is the percentile bootsrtap 1-alphav percent confidence interval
%          bootx is the boostrapped distribution is the estimator
%          p is the p value of the estimator to differ from mu
%
% GAR - University of Glasgow - Dec 2007
% added EST output - Nov 2008
% Cyril Pernet updated the function with rng and removed the bootstrap loop
%              for faster computation using randi , also allowed for matrix
%              input and added bootx output - Aug 2015
% ----------------------------------------------------------------------
% Copyright (C) RST Toolbox Team 2015

%% check arguments

rng('shuffle')

[p,N]=size(X); % number of estimates to compute
if p==1 % if X row vector, transpose to column
    X=X'; 
    [p,N]=size(X);
end

if nargin<2 || isempty(nboot)
  nboot=1000;
end

if nargin<3 || isempty(alphav)
    alphav=.05;
end
alphav = alphav / N;

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
        case ('median')
            q = [];
        otherwise
            error('the estimator name is not recognized, check arguments')
    end
end

%% compute the percentile CI

lo = round(nboot.*alphav./2);
hi = nboot - lo;

% cache
EST = NaN(1,N);
bootx = NaN(nboot,N);

for s = 1:N
    x = X(~isnan(X(:,s)),s);
    n = length(x);
    
    switch lower(est)
        case {'hd'} % Harrell-Davis estimator
            EST(s)=rst_hd(x,q);
            bootx(:,s)=rst_hd(x(randi(n,n,nboot)),q);
        case {'tm'} % trimmed mean
            EST(s)=rst_trimmean(x,q);
            bootx(:,s)=rst_trimmean(x(randi(n,n,nboot)),q);
        case {'trimmean'} % trimmed mean - Matlab function
            EST(s)=trimmean(x,q*2);
            bootx(:,s)=trimmean(x(randi(n,n,nboot)),q*2);
        otherwise
            EST(s)=median(x); % median
            bootx(:,s)=median(x(randi(n,n,nboot)));
    end
end

bootx = sort(bootx,1); % sort in ascending order
CI(1,:) = bootx(lo+1,:);
CI(2,:) = bootx(hi,:);

%% compute the null hypothesis
if nargout > 3
    if nargin < 6 || isempty(mu)
       mu=0; 
    end
   p = sum(bootx>mu,1)./nboot;
   p = 2.*min(p,1-p);
end
