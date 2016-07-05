function [diff,CI_boot,p,h,H]= rst_multicompare(Data,pairs,alphav,est,newfig,plothist)

% performs multiple comparisons between pairs based on a percentile bootstrap
% on differences - alphav is adjusted to control the type 1 error rate.
%
% FORMAT
% [diff,CI_boot,p,h]= rst_multicompare(Data,pairs,alphav,est,newfig,plothist)
%
% INPUT
% Data     = 2D matrix (n*p) of repeated measures
% pairs    = matrix 2*n for which pairs to test ([] means all)
% alphav    = 5% (default) or 10%
% est      = which estimator to use 'mean','median', or 'trimmean'
% newfig   = 0 if one doesn't want the function to create a new figure 
%            (ie allows to plot within an existing one) 
% plothist = 1 if one want to also plot the density histogram of each difference
%
% OUTPUT
% diff    = all pairwise differences between conditions
% CI_boot = percentile bootstrap confidence intervals of the differences
% p       = percentile bootstrap p values
% h       = significance of the p value (adjusted for multiple tests)
% H       = shows pairs and h
%
% Multiple comparisons are computed as described in Wilcox, R.R. (2005)
% Robust Estimation and Hypothesis Testing (2nd Ed). Elsevier, Academic
% Press, San Diego, CA, USA. If N<80 and less than 10 pairs are tested, the
% alphav value is set a priori based on simulation results. If N<80 and more
% then 10 pairs are tested, alphav is Bonferroni corrected. If N>80, alphav
% is adjusted using Hochberg step-up procedure. Overview of these
% procedures, and others, can be found at 
% http://en.wikipedia.org/wiki/Familywise_error_rate
%
% Cyril Pernet 16-08-2011
% ----------------------------------------------------------------------
% Copyright (C) RST Toolbox Team 2014


%% inputs

% set defaults
if nargin < 2
    pairs = [];
end

if nargin < 3
    alphav = 5/100;
end

if nargin < 4
    est = 'trimmean';
end

if nargin < 5
     newfig = 1;
end

if nargin < 6
     plothist = 0;
end

[n,p]=size(Data);
nboot = 1000;

% get pairs
% -----------
if isempty(pairs)
    all_pairs = nchoosek([1:p],2);
else
    all_pairs = pairs;
end
L = length(all_pairs);

% get alphav if n < 80
% -------------------
if n < 80 % small sample size, method SR
    if L>10
        alphav = alphav / L;
    else
        test(1) = (alphav ~= .05);
        test(2) = (alphav ~= .1);
        if sum(test) ==2
            error('multiple comparisons only works with apha 5% or 10%');
        else
            alphav = getalphavc(alphav,L);
        end
    end
end

low = round(nboot.*alphav./2);
high = nboot - low;

% compute bootstrap differences
% ------------------------------
boot_table = zeros(size(Data,1),nboot); B=1;
while B~=nboot+1
    tmp = ceil(rand(1,size(Data,1)).*size(Data,1));
    if length(unique(tmp)) ~= 1
        boot_table(:,B) = tmp;
        B=B+1;
    end
end

D = zeros(size(Data,1),L);
for i=1:L
    D(:,i) = Data(:,all_pairs(i,1)) - Data(:,all_pairs(i,2));
end 

% --------------------
if strcmp(est,'mean')
    diff = nanmean(D);
elseif strcmp(est,'median')
    diff = nanmedian(D);
elseif strcmp(est,'trimmean')
    diff = rst_trimmean(D); % 20% trim-mean
end

% --------------------

% bootstrap differences
avg_boot = zeros(nboot,L); 
if strcmp(est,'mean')
    for b=1:nboot
        avg_boot(b,:) = nanmean(D(boot_table(:,b),:));
    end
elseif strcmp(est,'median')
    for b=1:nboot
        avg_boot(b,:) = nanmedian(D(boot_table(:,b),:));
    end
elseif strcmp(est,'trimmean')
    for b=1:nboot
        avg_boot(b,:) = rst_trimmean(D(boot_table(:,b),:));
    end
end

avg_boot = sort(avg_boot);
% --------------------
CI_boot = [avg_boot(low,:);avg_boot(high,:)];
% --------------------

% get the stats
pl= sum((avg_boot>0)) / nboot;
% --------------------
p= 2*(min(pl,1-pl));
% --------------------
% check special case of 0 difference
test =sum(avg_boot);
for i=1:length(test)
    if test(i) == 0
        p(i) = 1;
    end
end

[sortedp,order]=sort(p','descend');
[reorder,inverseindex] = sort(order);

% return h
% ---------
h = zeros(L,1);  
if n < 80 % small sample size, method SR
    tmp = sortedp<alphav;
    h = tmp(inverseindex);
    
else % large sample size, Hochberg, 1988
    decreased_alphav = repmat(alphav,L,1)./(1:L)';
    test = sortedp < decreased_alphav;
    go = 1; index = 1; 
    while go == 1
        if test(index) == 1
            h(index:L) = 1;
            go = 0;
        elseif test(index) == 0
            index = index+1;
            if index == length(test)
                go = 0;
            end
        end
    end
    h = h(inverseindex);  
end

H = [all_pairs h];

% plot
if newfig == 1
    figure('Name', 'Pair-wise comparisons')
    set(gcf,'Color','w')
    bar([1:L],diff); hold on
    errorbar([1:L],diff,diff-CI_boot(1,:),CI_boot(2,:)-diff,'r','LineWidth',2);
    title(['Difference between pairs and ' num2str(100-alphav*100) '% CI'],'FontSize',16); grid on
    ylabel('Difference between condition','FontSize',14);
    for i=1:L
        xname{i} = [num2str(all_pairs(i,1)) '/' num2str(all_pairs(i,2)) ' '];
    end
    xlabel(sprintf('Pairs %s', cell2mat(xname)),'FontSize',14);
    set(gca,'FontSize',12);
    v=axis; ymin = v(3) + (v(4)-v(3))/10;
end

% additional plots?
if plothist == 1
    for i=1:L
        figure; rst_density_hist(D(:,i));
    end
end


end

function alphavc = getalphavc(alphav,L)

% table precomputed by Wilcox

if alphav == 5/100
    if L == 1; alphavc = .025;
    elseif L == 2; alphavc = .025;
    elseif L == 3; alphavc = .0169;
    elseif L == 4; alphavc = .0127;
    elseif L == 5; alphavc = .0102;
    elseif L == 6; alphavc = .00851;
    elseif L == 7; alphavc = .0073;
    elseif L == 8; alphavc = .00639;
    elseif L == 9; alphavc = .00568;
    elseif L == 10; alphavc = .00511;
    end
else
    if L == 1; alphavc = .005;
    elseif L == 2; alphavc = .005;
    elseif L == 3; alphavc = .00334;
    elseif L == 4; alphavc = .00251;
    elseif L == 5; alphavc = .00201;
    elseif L == 6; alphavc = .00167;
    elseif L == 7; alphavc = .00143;
    elseif L == 8; alphavc = .00126;
    elseif L == 9; alphavc = .00112;
    elseif L == 10; alphavc = .00101;
    end
end
end
