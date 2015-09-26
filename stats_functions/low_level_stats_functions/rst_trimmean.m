function [TM,CI] = rst_trimmean(varargin)

% This function computes the X percent trimmed mean of the data.
% It also returns the 95% confidence interval of the trimmed mean.
% 
% FORMAT
% [TM,CI] = rst_trimmean(data,percent)
%
% INPUT
% data is a vector or a matrix (also accept NaN)
% percent is the amount trimmed each side of data (default = 20%)
%
% OUTPUT
% TM is the trimmed mean value(s) of data taken column-wise
% CI is the 95% confidence interval
%
% the trimmed mean is the mean of the data excluding the highest and lowest
% K data values (K=round(N*percent)) and N is the number of values
% in data) - see e.g. Wilcox, Rand R. "Introduction to Robust Estimation
% and Hypothesis Testing." New York: Academic Press. 2005.
%
% Cyril Pernet v1 - April 2012
% -------------------------------------------------------------------------
% Copyright (C) RST Toolbox Team 2015

%% check arguments and data

data = varargin{1};
percent = 20/100; % default

if nargin>3
    error('too many arguments')
elseif nargin == 2
    percent = varargin{2};
    if ~isnumeric(percent)
            error('percent of trimming must be a numeric')
    elseif percent > 1
        percent = percent / 100;
    end
end

[n,p]=size(data);
if n== 1 && p>2 
    data = data'; % transpose row vector in column vector if needed
    [n,p]=size(data);
end


%% down to business
for c=1:p
    
    tmp = data(:,c);
    tmp(isnan(tmp)) = [];
    N = size(tmp,1);
    K=round(N*percent);
    if K == 0
        error('not enough data points to trim')
    end
    X = sort(tmp,1);
    TM(c) = mean(X((K+1):(N-K),:),1);
    
    %% add robust CI
    if nargout > 1
        [tv,g] = tvar(tmp,percent*100); % trimmed squared standard error + g trimmed elements
        se = sqrt(tv); % trimmed standard error
        df = N - 2.*g - 1; % n-2g = number of observations left after trimming
        CI(1,c) = TM(c)+tinv(percent./2,df).*se;
        CI(2,c) = TM(c)-tinv(percent./2,df).*se;
    end 
end

end

function [tv,g]=tvar(x,percent)

% function [tv]=tvar(x,percent)
% The trimmed variance is calculated from the winsorized variance, vw,
% according to the formula vw/(k*length(x)), where k=(1-2*percent/100)^2.
% Original code provided by Prof. Patrick J. Bennett, McMaster University
% See Rand R. Wilcox (2001), Fundamentals of Modern Statisical Methods, page 164.
% See also Rand R. Wilcox (2005), p.61-63
% Edit input checks: GA Rousselet - University of Glasgow - Nov 2008
% Merge function tvar, winvar, winsample - C Pernet June 2010

g=floor((percent/100)*size(x,1));
xsort=sort(x,1);
loval=xsort(g+1,:);
hival=xsort(size(x,1)-g,:);
for i=1:size(x,2)
    x(find(x(i,:)<=loval(i)),i) = loval(i); % instead of trimming like in TM(c) we substutute values
    x(find(x(i,:)>=hival(i)),i) = hival(i);
    wv(i)=var(x(:,i));
end
k=(1-2*percent/100)^2;
tv=wv/(k*length(x));

end



