function [diff,CI,p,alphav,h]= rst_multicompare(Data,Pairs,varargin)

% Performs multiple comparisons between pairs (repeated measures) based on 
% a percentile bootstrap on differences - the alpha value is adjusted to 
% control the type 1 error rate. This is not the same as post-hoc testing
% in ANOVA in which significant effect across all conditions was found in
% a 1st place. The control is across tested difference only. Multiple 
% comparisons are computed as described in Wilcox, R.R. (2012). Robust 
% Estimation and Hypothesis Testing (3rd Ed). Elsevier, Academic Press, 
% San Diego, CA, USA. p398-401 (section 8.3.2).
%
% FORMAT
% [diff,CI,p,alphav,h]= rst_multicompare(Data,Pairs,'alphav',0.05,'estimator','median','newfig','yes')
%
% INPUT
% Data      = 2D matrix (n*p) of repeated measures
% pairs     = matrix 2*n for which pairs to test ([] means all)
% alphav    = 5% (default) - for n<80 and the number of pairs is more than
%             10, only 5% or 10% are allowed
% estimator = which estimator to use 'mean','median' (Harell-Davis) or 'trimmed mean'
% newfig    = 'yes' or 'no' allows to plot within an existing one 
%
% OUTPUT
% diff     = all pairwise differences between conditions
% CI       = percentile bootstrap confidence intervals of the differences
% p        = the percentile bootstrap p values
% alphav   = the adjusted alpha value used for significance
% h        = significance of the p value (adjusted for multiple tests)
%
% usage
% Data = randn(30,6); Data(:,5) = Data(:,5)+1; % no differences but column 5
% Pairs = [1 4; 2 5; 3 6];
% [diff,CI,p,alphav,h]= rst_multicompare(Data,Pairs,'alphav',0.05,'estimator','trimmed mean','newfig','yes')
% 
% The tested hypothesis is that each diffence scores Dij = 0. This is
% performed by bootstrapping the matrix of differences, computing N estimators
% and testing how many times it is bigger than 0. FWE rate is controlled 
% using the SR method if N<80 and less than 10 pairs are tested. If N<80 
% and more than 10 pairs are tested, alphav adjusted for the number of tests.
% If N>80, alphav is adjusted using Hochberg step-up procedure. 
% Overview of these procedures, and others, can be found at 
% <http://en.wikipedia.org/wiki/Familywise_error_rate>
%
% Cyril Pernet
% ----------------------------------------------------------------------
% Copyright (C) RST Toolbox Team 2017

% extra input for quality check
% 'plothist' if one wants to plot the density histograms of each difference


%% set defaults

% _hard coded_ 
nboot = 1000; % how many bootstraps are used

% _soft coded_
alphav    = 5/100; % non adjusted alpha value
estimator ='median'; % by default use a robust estimator
newfig    = 'yes'; % start a new figure
plothist  = 'no'; % do not create additional histograms of differences

% _check inputs_ 
if size(Data,2) < 2
    error('The Data matrix must be at least 2 column to make a pair')
else
  n = size(Data,1);    
end

% compute pairs
if nargin ==1 
    all_pairs = nchoosek([1:size(Data,2)],2);
else
    if size(Pairs,2) ~= 2
        error('Pairs must be a n*2 matrix')
    end
    all_pairs = Pairs; clear Pairs;
end
L = length(all_pairs);

% check options
for in=1:length(varargin)
        if strcmpi(varargin(in),'alphav')
            alphav = cell2mat(varargin(in+1));
            if alphav > 1; alphav = alphav / 100; end
        elseif strcmpi(varargin(in),'estimator')
            if ~strcmpi(varargin(in+1),'median') && ...
                    ~strcmpi(varargin(in+1),'mean') && ...
                    ~strcmpi(varargin(in+1),'trimmed mean')
                error(['estimator ' cell2mat(varargin(in+1)) ' is not recognized'])
            else
                estimator = cell2mat(varargin(in+1));
            end
        elseif strcmpi(varargin(in),'newfig')
            newfig = cell2mat(varargin(in+1));
        elseif strcmpi(varargin(in),'plothist')
            plothist = cell2mat(varargin(in+1));
        end
end

% non-adjusted confidence intervals
low = round(nboot.*alphav./2);
high = nboot - low;
orig_alpha = alphav;

% get the adjusted alpha value
if n < 80 && L<=10 % small sample size, method SR
        if sum([alphav ~= .05 alphav ~= .1]) == 2
            error('multiple comparisons with N<80 and less than 10 pairs, only works with apha 5% or 10%');
        else
            tmp = nan(1,L);
            for t=1:L
                tmp(t) = getalphavc(alphav,t); % see sub-function getalphav 
            end
            alphav = tmp; clear tmp
        end
else % small samples and more than 10 tests or large sample: use Hochberg method
    alphav = repmat(alphav,1,L) ./ [1:L];
end

%% start the bootstrap

% _compute the differences_  
D = zeros(size(Data,1),L);
for i=1:L
    D(:,i) = Data(:,all_pairs(i,1)) - Data(:,all_pairs(i,2));
end 

if strcmpi(estimator,'mean')
    diff = nanmean(D);
elseif strcmpi(estimator,'median')
    diff = rst_hd(D,0.5); % use the Harell David estimator
elseif strcmpi(estimator,'trimmed mean')
    diff = rst_trimmean(D); % 20% trimmed mean
end

% _compute bootstrap differences_  

% make a boot table with at least 2 different observations  
boot_table = randi(n,[n nboot*2]);
for B =1:nboot*2
    if length(unique(boot_table(:,B))) < 3 
        boot_table(:,B) = NaN;
    end
end
boot_table(:,isnan(sum(boot_table,1))) = [];
boot_table = boot_table(:,1:nboot);

% _do all the bootstraps at once_  
avg_boot = zeros(nboot,L); 
for d=1:L % for each pair
    tmp = D(:,d);
    if strcmpi(estimator,'mean') % compute the nboot estimators
        avg_boot(:,d) = nanmean(tmp(boot_table));
    elseif strcmpi(estimator,'median')
        avg_boot(:,d) = rst_hd(tmp(boot_table),0.5);
    elseif strcmpi(estimator,'trimmed mean')
        avg_boot(:,d) = rst_trimmean(tmp(boot_table));
    end
end
avg_boot = sort(avg_boot,1);
CI       = [avg_boot(low,:);avg_boot(high,:)]; % get the confidence intervals

%% get the p values and h

pl= sum((avg_boot>0)) / nboot;
p= 2*(min(pl,1-pl)); % generalized p-value
p(p==0) = 1/nboot; % check p-values of 0, change to precision
p(sum(avg_boot,1)==0) = 1; % check special case of 0 difference for all pairs, change to 1

% return h based on adjusted alpha-value
h = zeros(L,1);  
[sortedp,order]= sort(p','descend');
[~,inverseindex] = sort(order);
for t=1:L
    test = sortedp(t) <= alphav(t);
    if test == 1
        h(t:L) = 1; break
    end
end
h = h(inverseindex);

%% figures

if strcmpi(newfig,'yes')
    color_scheme = cubehelixmap('semi_continuous',L+10);
    figure('Name', 'Pair-wise comparisons');set(gcf,'Color','w')
    for i=1:L
        rectangle('Position',[i-0.2,CI(1,i),0.4,CI(2,i)-CI(1,i)],'Curvature',[0.4 0.4],'LineWidth',2,...
            'FaceColor',color_scheme(i+8,:),'EdgeColor',[0.35 0.35 0.35]); hold on;
        plot([i-0.2 i+0.2],[diff(i) diff(i)],'LineWidth',3,'Color',[0.35 0.35 0.35]);
        tmp = repmat(D(:,i),1,2); tmp(1:2:end,1) = NaN; tmp(2:2:end,2) = NaN; 
        scatter(repmat(i-0.05,n,1),tmp(:,1),50,'k','Filled');
        scatter(repmat(i+0.05,n,1),tmp(:,2),50,'k','Filled');
    end
    title(sprintf('%s difference between \n pairs with %g %% CI',estimator,100-orig_alpha*100),'FontSize',16); 
    ylabel('Difference between pairs');
    for i=1:L
        xname{i} = [num2str(all_pairs(i,1)) '/' num2str(all_pairs(i,2)) ' '];
    end
    xlabel(sprintf('Pairs %s', cell2mat(xname))); box on; grid on
    set(gca,'FontSize',12); v=axis; axis([v(1) v(2) v(3)-0.5 v(4)+0.5]);
end

if strcmpi(plothist,'yes')
    for i=1:L
        figure; rst_density_hist(avg_boot(:,i)); 
        title(['bootstrap difference pair ' num2str(i)]);
    end
end


end

function alphavc = getalphavc(alphav,L)

% table precomputed by Wilcox p347

if alphav == 5/100
    if     L == 1;  alphavc = .025;
    elseif L == 2;  alphavc = .025;
    elseif L == 3;  alphavc = .0169;
    elseif L == 4;  alphavc = .0127;
    elseif L == 5;  alphavc = .0102;
    elseif L == 6;  alphavc = .00851;
    elseif L == 7;  alphavc = .0073;
    elseif L == 8;  alphavc = .00639;
    elseif L == 9;  alphavc = .00568;
    elseif L == 10; alphavc = .00511;
    end
else
    if     L == 1;  alphavc = .005;
    elseif L == 2;  alphavc = .005;
    elseif L == 3;  alphavc = .00334;
    elseif L == 4;  alphavc = .00251;
    elseif L == 5;  alphavc = .00201;
    elseif L == 6;  alphavc = .00167;
    elseif L == 7;  alphavc = .00143;
    elseif L == 8;  alphavc = .00126;
    elseif L == 9;  alphavc = .00112;
    elseif L == 10; alphavc = .00101;
    end
end
end
