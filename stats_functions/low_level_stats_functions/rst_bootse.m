function bse = rst_bootse(x,nboot,est,q)

% See Wilcox p.44-46
%
% GAR, University of Glasgow, Dec 2007

rand('state',sum(100*clock));

if nargin<2;nboot=1000;est='median';end
if nargin<4;q=.5;end

n=length(x); % number of subjects
boot = zeros(1,nboot);

switch lower(est)
    case {'hd'} % Harrell-Davis estimator
        for kk=1:nboot % do bootstrap
            eval(['boot(kk)=',est,'(randsample(x,n,true),q);'])
        end
    otherwise
        for kk=1:nboot % do bootstrap
            eval(['boot(kk)=',est,'(randsample(x,n,true));'])
        end
end

% We estimate the standard error of the est statistic by the standard
% deviation of the bootstrap replications
% Efron & Tibshinari 1993, chapter 6
% Wilcox 2005, p.44-45
bse = std(boot,0); % normalize by (n-1)

return
