function [HDQ,CIQ] = rst_hd(varargin)

% Compute the Harrell-Davis estimate of the qth decile
% The vector x contains the data, and the desired decile is q
%
% FORMAT [HDQ,CIQ] = rst_hd(x,q,nboot)
%
% INPUT  x is a vector or a matrix, in this case HD estimates are computed column-wise
%          note that matrices can contain NaN - and this is accounted for
%        q is the decile to compute (default is the median q=0.5)
%        nboot the number of bootstrap to use to compute the 95% CI (default = 100)
%
% OUTPUT HDQ is/are the decile(s)
%        CIQ is/are the 95% confidence interval
%
% FRANK E. HARRELL and C. E. DAVIS (1982).
% A new distribution-free quantile estimator
% Biometrika 69 (3): 635-640. doi: 10.1093/biomet/69.3.635
%
% Original R code by Rand Wilcox, p.71, 139 and p.130-134
% GAR, University of Glasgow, Dec 2007 / Sep 2009
% CP, The University of Edinburgh - update of the code combining estimate
% and CI + takes matrices as input August 2014
% ----------------------------------------------------------------------
% Copyright (C) RST Toolbox Team 2015

%% inputs

% defaults
X = varargin{1};
[p,N]=size(X); % number of estimates to compute
if p==1 % if X row vector, transpose to column
    X=X'; 
    [p,N]=size(X);
end
q=.5; % median
nboot = 1000; % use for 95% CI

if nargout == 1 && nargin > 2
    disp('too many arguments in, only the 2 first ones will be used')
elseif nargout == 1 && nargin == 2
    q = varargin{2};
elseif nargout == 2
    if p<=10
        error('confidence intervals of the hd estimates of the deciles cannot be computed for less than 11 observations')
    end
    rng shuffle
    if nargin >1; q = varargin{2}; end
    if nargin >2; nboot = varargin{3}; end
end

%% compute
for i=1:N
    x = X(~isnan(X(:,i)),i);
    HDQ(i) = get_HD(x,q);
    
    if nargout == 2
        % The constant c was determined so that the
        % probability coverage of the confidence interval is
        % approximately 95% when sampling from normal and
        % non-normal distributions
        n=length(x);
        c = 1.96 + .5064.* (n.^-.25);
        
        if n<=10
            error(sprintf('confidence intervals of the hd estimates of the deciles cannot be computed \n for less than 11 observations, error column %g',i))
        elseif n <=21 &&  q<=.2 || n <=21 && q>=.8
            c = -6.23./n+5.01;
        elseif n<=40 && q<=.1 || n<=40 && q>=.9
            c = 36.2./n+1.31;
        end
        
        % do bootstrap
        for kk=1:nboot
            boot(kk)=get_HD(randsample(x,n,true),q);
        end
        bse = std(boot,0); % normalize by (n-1)
        CIQ(1,i) = HDQ(i)-c.*bse;
        CIQ(2,i) = HDQ(i)+c.*bse;
    end
    clear x n c boot bse
end

end

function HD = get_HD(x,q)
% that's the Harrell-Davis estimator
n=length(x);
m1=(n+1).*q;
m2=(n+1).*(1-q);
vec=1:n;
w=betacdf(vec./n,m1,m2)-betacdf((vec-1)./n,m1,m2);
y=sort(x);
HD=sum(w(:).*y(:));
clear vec w y
end


