function [tmean,tmCIlo,tmCIhi,df,se,p,Ty] = rst_trimci(x,percent,alpha,out,mu)

% function [tm,tmCIlo,tmCIhi,g,se,p,Ty] = trimci(x,percent,alpha,out,mu)
% Returns the 1-alpha confidence interval of the trimmed mean.
% 
% INPUTS:
%         x - vector
%         percent - percentage of trimming, default = 20
%         alpha - default = .05
%         out - output results in command window, default 1
%         
% OUTPUTS:
%         tm - trimmed mean
%         tmCIlo - lower bound of the CI
%         tmCIhi - upper bound of the CI
%         g - number of trimmed data points on each side of the distribution
%         se - estimate of the standard error of the trimmed mean
%         p - probability of obtaining a t-value whose absolute value
%             is greater than Ty computed under the null hypothesis, by default mu=0
%         Ty - t value for trimmed sample
% 
% GAR, University of Glasgow, 23 Nov 2008
% See Rand R. Wilcox (2005), p.113-117
%
% See also WINVAR, TVAR

if nargin < 2
    percent = 20;
    if nargin < 3
        alpha = .05;
        if nargin < 4
            out = 1;
            if nargin < 5
                mu = 0;
            end
        end
    end
end

if percent >= 100 || percent < 0
    error('trimci:InvalidPercent', 'PERCENT must be between 0 and 100.');
end

% make sure that x is a vector
sz = size(x);
if sz > 2 | min(sz) > 1    
  error('trimci requires x to be a vector, not a matrix.');
end    

tmean = rst_trimmean(x,percent); % compute trimmed mean
[tv,g] = rst_winvar(x,percent); % trimmed squared standard error + g trimmed elements
se = sqrt(tv); % trimmed standard error
df = length(x) - 2.*g - 1; % n-2g = number of observations left after trimming

% calculate CI for trimmed mean
tmCIlo = tmean+tinv(alpha./2,df).*se; % 
tmCIhi = tmean-tinv(alpha./2,df).*se; % 

if out==1 % optional command window output
fprintf('\n-----------------------------------------------------\n')
fprintf('%g percent trimmed mean %g percent confidence interval:\n',percent,(100.*(1-alpha)))
fprintf('trimmed mean: %5.2f < %5.2f < %5.2f\n',tmCIlo,tmean,tmCIhi);
fprintf('-----------------------------------------------------\n')
end

if nargout > 5
    Ty = (tmean-mu)./se;
    p=2*(1-tcdf(abs(Ty),df));
    h=abs(Ty)>tinv(1-alpha,df);
end


