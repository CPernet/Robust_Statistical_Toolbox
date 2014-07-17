function [regout,CI,pboot,p] = rst_pbCIreg(a,b,nboot,alpha,mu)

% function [regout,CI,pboot,p] = pbCIreg(a,b,nboot,alpha,mu)
%
% Compute a percentile bootstrap confidence interval
% of the linear regression outputs using OLS
%
%   INPUTS:
%           a (vector of predictor)
%           b (vector of response)
%           nboot (resamples, default 600)
%           alpha (default 5% - i.e. type .05)
%           mu (optional null hypothesis, necessary to get a p-value in
%           one-sample hypothesis testing procedure, default 0)
%
% GAR - University of Glasgow - Aug 2008

rand('state',sum(100*clock));

if nargin<3 || isempty(nboot)
    nboot=600;
end
if nargin<4 || isempty(alpha)
    alpha=.05;
end

if length(a) ~= length(b)
    error('error, input vectors must have same size')
end

% check vector orientation
ss=size(a);
if ss(1)<ss(2)
    a=a';
end

ss=size(b);
if ss(1)<ss(2)
    b=b';
end

% create data matrix
n = length(a);
% a = cat(2,ones(n,1),a);

% do regression & get outputs
[regout,bint,r,rint,stats]=regress(a,[ones(length(b),1) b]);
p=stats(3);
% [B,BINT,R,RINT,STATS]=regress(b,a);
% slope=B;
% r2=STATS(1);
% p=STATS(3);

% set confidence interval boundaries
if nboot~=600
    fprintf('\nBetter confidence intervals can be obtained by using 600 random samples')
    lo = round(nboot.*alpha./2);
    hi = nboot - lo;
elseif n<40
    lo = 6;
    hi = 593;
elseif n>=40 && n<80
    lo = 7;
    hi = 592;
elseif n>=80 && n<180
    lo = 10;
    hi = 589;
elseif n>=180 && n<250
    lo = 13;
    hi = 586;
elseif n>=250
    lo = 15;
    hi = 584;
end

% do bootstrap
bootreg=zeros(nboot,size(b,2)+1);
for kk = 1:nboot % sample with replacement loop
    bootsample=randsample(1:n,n,'true');
    %     bootreg(kk)=regress(b(bootsample),a(bootsample,:));
    bootreg(kk,:)=regress(a(bootsample),[ones(length(b),1) b(bootsample,:)]);
end

diffsort = sort(bootreg); % sort in ascending order

for sub=1:length(regout)
    CI(sub,1) = diffsort(lo+1,sub);
    CI(sub,2) = diffsort(hi,sub);
end

% if nargout > 2
if nargin < 4 || isempty(mu)
    mu=0;
end

for sub=2:length(regout)
    pboot(sub) = sum(bootreg(:,sub)>mu)./nboot;
    pboot(sub) = 2.*min(pboot(sub),1-pboot(sub));
end
% end

return