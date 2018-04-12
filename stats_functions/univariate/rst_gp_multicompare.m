function [diff,CI,p,alphav,h]=rst_gp_multicompare(gp1,gp2,varargin)

% Test for differences between groups with an adjustment for multiple 
% comparisons - This is not the same as post-hoc testing in ANOVA in which 
% a significant interaction was found in a 1st place. The control is across
% tested difference only. The control is across tested difference only. 
% Multiple comparisons are computed as described in Wilcox, R.R. (2012). 
% Robust Estimation and Hypothesis Testing (3rd Ed). Elsevier, Academic 
% Press, San Diego, CA, USA. p329-330 (section 7.4.7).
%
% FORMAT [diff,CI,p,alphav,h]=rst_gp_multicompare(gp1,gp2,'alphav',0.05,'estimator','median','newfig','yes')
%
% INPUT gp1 and gp2 are 2 matrices n*p and m*p 
%                   with p the conditions to compare column-wise 
%       alphav    = 5% (default) - for n<80 and the number of pairs is more than
%                   10, only 5% or 10% are allowed
%       estimator = which estimator to use 'mean','median' (Harell-Davis) or 'trimmed mean'
%       newfig    = 'yes' or 'no' allows to plot within an existing one 
%
% OUTPUT difference is the estimator difference
%        CI the adjusted 95% CI
%        p        = the percentile bootstrap p values
%        alphav   = the adjusted alpha value used for significance
%        note there is no p values when comparing means
%        h  = significance of the test
%
% usage
% gp1 = randn(30,3); gp2 = randn(30,3); gp2(:,2) = gp2(:,2)+1; % no differences but gp2
% [diff,CI,p,alphav,h]= rst_gp_multicompare(gp1,gp2,'alphav',0.05,'estimator','trimmed mean','newfig','yes')
% 
% The tested hypothesis is that each group diffence Dij = 0. 
% When testing mean differences, adjustement is based on the maximum 
% differences over all bootstraps across conditions using a bootstrap t.
% This is limited to N<80. Note that there is no p values associated.
% For other estimators, FWE rate is controlled using the SR method 
% if N<80 and less than 10 pairs are tested. If N<80 and more than 10 pairs
% are tested, alphav adjusted for the number of tests. If N>80, alphav is 
% adjusted using Hochberg step-up procedure. Overview of these procedures, 
% and others, can be found at <http://en.wikipedia.org/wiki/Familywise_error_rate>
%
% Cyril Pernet 
% ----------------------------------------------------------------------
% Copyright (C) RST Toolbox Team 2017

% extra input for quality check
% 'plothist' if one wants to plot the density histograms of each difference

% _hard coded_ 
nboot = 1000; % how many bootstraps are used

% _soft coded_
alphav   = 5/100; % non adjusted alpha value
estimator ='median'; % by default use a robust estimator
newfig    = 'yes'; % start a new figure
plothist  = 'no'; % do not create additional histograms of differences

% _check inputs_ 
if any(size(gp1,2) ~= size(gp2,2))
    error('input data must have the same number of variables')
else
  [n,p] = size(gp1);    
end

% check options
for in=1:length(varargin)
    if strcmpi(varargin(in),'alphav')
        alphav = cell2mat(varargin(in+1));
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
orig_alpha = alphav;
low = round(nboot.*orig_alpha./2);
high = nboot - low;

%% multiple comparisons between 2 group means
if strcmpi(estimator,'mean')
    
    if length(gp1)>80 && length(gp2) > 80
        disp('adjustment needed for large N not implemented (see Wilcox 2nd edition, p314'); return
    end
    
    boot_gp1 = zeros(size(gp1,1),size(gp1,2),nboot); B=1;
    while B~=nboot+1
        tmp1 = randi(size(gp1,1),1,size(gp1,1))';
        if length(unique(tmp1)) ~= 1
            boot_gp1(:,:,B) = gp1(tmp1,:);
            B=B+1;
        end
    end
    
    boot_gp2 = zeros(size(gp2,1),size(gp2,2),nboot); B=1;
    while B~=nboot+1
        tmp2 = randi(size(gp2,1),1,size(gp2,1))';
        if length(unique(tmp2)) ~= 1
            boot_gp2(:,:,B) = gp2(tmp2,:);
            B=B+1;
        end
    end
    
    means_gp1 = squeeze(nanmean(boot_gp1,1)); % means over resampled subjects
    means_gp2 = squeeze(nanmean(boot_gp2,1));
    boot_means1 = nanmean(means_gp1,2); % means over all bootstrap
    boot_means2 = nanmean(means_gp2,2);
    Tau1 = (1/(nboot-1))*nansum((means_gp1-boot_means1(:,ones(1,nboot))).^2,2); % standard error of each condition
    Tau2 = (1/(nboot-1))*nansum((means_gp2-boot_means2(:,ones(1,nboot))).^2,2);
    H = abs(means_gp1 - means_gp2 - repmat((boot_means1-boot_means2),1,nboot)) ...
        ./ repmat((sqrt(Tau1+Tau2)),1,nboot); % difference between conditions
    Hmax = sort(max(H)); % maximum differences over all bootstraps across conditions
    
    u = round((1-alphav)*nboot);
    diff = nanmean(gp1,1)' - nanmean(gp2,1)';
    bound = Hmax(u).*(sqrt(Tau1+Tau2)); % because of Hmax, this adjusts for multiple testing
    CI = NaN(length(diff),2);
    CI(:,1) = diff - bound;
    CI(:,2) = diff + bound;
    h = single(CI(:,1) > 0) + single(CI(:,2) < 0);
    p = NaN; CI = CI';
end


%% multiple comparisons between 2 group medians or trimmed means
if ~strcmpi(estimator,'mean')
    
    % adjusted alpha value
    if n < 80 && p<=10 % small sample size, method SR
        if sum([alphav ~= .05 alphav ~= .1]) == 2
            error('multiple comparisons with N<80 and less than 10 pairs, only works with apha 5% or 10%');
        else
            tmp = nan(1,p);
            for t=1:p
                tmp(t) = getalphavc(alphav,t); % see sub-function getalphav
            end
            alphav = tmp; clear tmp
        end
    else % small samples and more than 10 tests or large sample: use Hochberg method
        alphav = repmat(alphav,1,p) ./ [1:p];
    end
    
    boot_gp1 = zeros(size(gp1,1),size(gp1,2),nboot); B=1;
    while B~=nboot+1
        tmp1 = randi(size(gp1,1),1,size(gp1,1))';
        if length(unique(tmp1)) ~= 1
            boot_gp1(:,:,B) = gp1(tmp1,:);
            B=B+1;
        end
    end
    
    boot_gp2 = zeros(size(gp2,1),size(gp2,2),nboot); B=1;
    while B~=nboot+1
        tmp2 = randi(size(gp2,1),1,size(gp2,1))';
        if length(unique(tmp2)) ~= 1
            boot_gp2(:,:,B) = gp2(tmp2,:);
            B=B+1;
        end
    end
    
    % compute the difference between estimators
    D = NaN(p,nboot);
    CI = NaN(2,p);
    for t=1:p
        if strcmpi(estimator,'median')
            diff(t) = rst_hd(gp1(:,t),0.5) - rst_hd(gp2(:,t),0.5);
            boot_est1 = rst_hd(squeeze(boot_gp1(:,t,:)),0.5);  % use the Harell David estimator
            boot_est2 = rst_hd(squeeze(boot_gp2(:,t,:)),0.5);
        elseif strcmpi(estimator,'trimmed mean')
            diff(t) = rst_trimmean(gp1(:,t)) - rst_trimmean(gp2(:,t));
            boot_est1 = rst_trimmean(squeeze(boot_gp1(:,t,:))); % 20% trimmed mean
            boot_est2 = rst_trimmean(squeeze(boot_gp2(:,t,:)));
        end
        D(t,:) = sort(boot_est1 - boot_est2);
        CI(:,t) = [D(t,low) D(t,high)];
    end
    
    % get the p values and h
    pl= sum((D>0),2) / nboot;
    p= 2*(min(pl,1-pl)); % generalized p-value
    p(p==0) = 1/nboot; % check p-values of 0, change to precision
    p(sum(D,2)==0) = 1; % check special case of 0 difference for all pairs, change to 1
    
    % return h based on adjusted alpha-value
    L = size(gp1,2); h = zeros(L,1);
    [sortedp,order]= sort(p','descend');
    sortedalphav = alphav(order);
    [~,inverseindex] = sort(order);
    
    for t=1:L
        test = sortedp(t) <= sortedalphav(t);
        if test == 1
            h(t:L) = 1; break
        end
    end
    h = h(inverseindex);
end


%% figure
if strcmpi(newfig,'yes')
    figure('Name', 'group comparisons'); L = size(gp1,2); 
    set(gcf,'Color','w'); color_scheme = cubehelixmap('semi_continuous',L+10);
    
    if size(gp1) == size(gp2)
        D = gp1 - gp2; plot_D = 1;
    else
    plot_D = 0;
    end
    for i=1:L
        rectangle('Position',[i-0.2,CI(1,i),0.4,CI(2,i)-CI(1,i)],'Curvature',[0.4 0.4],'LineWidth',2,...
            'FaceColor',color_scheme(i+8,:),'EdgeColor',[0.35 0.35 0.35]); hold on;
        plot([i-0.2 i+0.2],[diff(i) diff(i)],'LineWidth',3,'Color',[0.35 0.35 0.35]);
        if plot_D == 1
            tmp = repmat(D(:,i),1,2); tmp(1:2:end,1) = NaN; tmp(2:2:end,2) = NaN; 
        scatter(repmat(i-0.05,n,1),tmp(:,1),50,'k','Filled');
        scatter(repmat(i+0.05,n,1),tmp(:,2),50,'k','Filled');
        end
    end
    if strcmpi(estimator,'mean')
        title('Mean differences between groups \n with adjusted CI','FontSize',16); 
    else
        title(sprintf('%s differences between \n groups with %g  CI',estimator,100-orig_alpha*100),'FontSize',16); 
    end

    ylabel('Difference between groups');box on; grid on
    xlabel(sprintf('Test %g ',[1:L])); 
    set(gca,'FontSize',12); v=axis; axis([v(1) v(2) v(3)-0.5 v(4)+0.5]);
end

end

function alphavc = getalphavc(alphav,L)

% table precomputed by Wilcox p330

if alphav == 5/100
    if     L == 1;  alphavc = .05;
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
    if     L == 1;  alphavc = .01;
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


