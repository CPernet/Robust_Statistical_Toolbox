function [bse,CI] = rst_bootse(varargin)

% Computes the bootstrap standard error
% 
% FORMAT [bse,CI] = rst_bootse(X,estimator,nboot,q)
%
% INPUT  X is a vector or a matrix, in this case HD estimates are computed column-wise
%          note that matrices can contain NaN - and this is accounted for
%        estimator should be 'mean', 'trimmed mean', 'median', 'Harrell-Davis'
%        nboot is the number of bootstrap samples (default = 1000)
%        q is either the 'Harrell-Davis' quantile (default is q=0.5, estimate 
%        of the median) or the percentage of trimming for 'trimmed mean'
%        (default 20%)
%
% OUTPUT bse is the standard error (ie the standard deviation of the
%                                   bootstrapped estimator)
%        CI is the 95% Confidence Interval of the estimator computed for
%        each column (non adjusted for multiple testing)
%
% See Wilcox p.44-46
% GAR, University of Glasgow, Dec 2007
% CP, The University of Edinburgh, Reformated for many estimators, 
% removed the bootstrap loop for a resampling table/vectorization and
% added the 95% CI output - 08/2014
% -------------------------------------------------------------------------
% Copyright (C) RST Toolbox Team 2014

%% inputs

% defaults
X = varargin{1};
est = varargin{2};
[p,N] = size(X);
if p == 1 && N > 2
    X = X';
    [p,N] = size(X);
end
nboot = 1000;

% check estimator
EST{1} = 'mean'; EST{2} = 'median'; EST{3} = 'trimmed mean'; EST{4} = 'harrell-davis';
if ~ischar(est)
    error('estimator argument must be a character')
elseif sum(strcmp(lower(est),EST)) ==0
    error('unknown estimator name')
end

if nargin == 2 || nargin == 3
    if strcmp(lower(est),'harrell-davis')
        q = 0.5; % median
    elseif strcmp(lower(est),'trimmed mean')
        q = 0.2; % 20% trimming
    end
elseif nargin == 4 && strcmp(lower(est),'harrell-davis') || nargin == 4 && strcmp(lower(est),'trimmed mean')
    q = varargin{4};
end

bse = NaN(1,N);
if nargout == 2
    low = round((0.05*nboot)/2);
    high = nboot - low;
    CI = NaN(2,N);
end

%% compute column wise (allows taking care of NaNs and vectorizing the bootstrap)
for i=1:N
    x = X(~isnan(X(:,i)),i);
    n=length(x); % number of subjects
    boot = zeros(1,nboot);
    % make a resampling table with at least 2 different values
    table = NaN(n,nboot); b = 1;
    while b<nboot+1
        tmp = randi(n,[n,1]);
        if length(unique(tmp))>=2
            table(:,b) = tmp;
            b = b +1;
        end
    end
    
    % simply compute
    if strcmp(lower(est),'harrell-davis')
        boot= rst_hd(x(table),q);
    elseif strcmp(lower(est),'trimmed mean')
        boot= rst_trimmean(x(table),q);
    elseif strcmp(lower(est),'median')
        boot= median(x(table));
    elseif strcmp(lower(est),'mean')
            boot= mean(x(table));
    end  
    
    % We estimate the standard error of the est statistic by the standard
    % deviation of the bootstrap replications - Efron & Tibshinari 1993, chapter 6
    % Wilcox 2005, p.44-45
    bse(i) = std(boot,0); % normalize by (n-1)
    if nargout == 2
        boot = sort(boot);
        CI(1,i) = boot(low);
        CI(2,i) = boot(high);
    end
    clear x n boot table
end
