function [h,CI,p] = rst_pttest(varargin)

% performs a robust paired t test, ie tests if an estimator of the
% difference between pairs differs from zeros.
%
% FORMAT
% [h,CI] = rst_pttest(X,Y)
% [h,CI,p] = rst_pttest(X,Y,estimator,plot_option,alpha,nboot)
%
% INPUTS
% - X,Y are either vectors or a matrices, in the latter case pairs are
% taken column wise and the resampling is identical in each column
% - estimator = 'mean' --> performs a bootstrap-t on the mean
% - estimator = 'median' --> performs a percentile bootstrap on the median
% - estimator = 'trimmean' --> performs a percentile bootstrap on the 20% trimmean (default)
% - plot_option indicates to plot the result (1 - default) or not (0)
% - alpha is the alpha level, default = 5%
% - nboot is the number of bootstrap, default = 600
%
% OUTPUTS
% [h,CI] = rst_pttest(X,Y,'mean')
% h is the hypothesis that the mean difference is 0 if h = 0 then this is not
% significantly different from 0, i.e. pairs are not different 
% note that if X and Y are matrices, a Bonferonni correction is applied for the number
% of tests (columns)
% CI are the alpha percent (Bonferoni corrected) confidence intervals used to get h
%
% [h,CI,p] = rst_pttest(data,'median');
% [h,CI,p] = rst_ptest(data,'trimmean');
% h is the hypothesis that the median difference or 20% trimmean difference is 0 
% if h = 0 then this is not significantly different from 0, i.e. pairs do not differ 
% note that if X and Y are matrices, a Bonferonni correction is applied for the number
% of tests (columns) 
% p are the p values
% CI are the alpha percent (Bonferoni corrected) confidence intervals
%
% Cyril Pernet - v1 - April 2012
% -------------------------------

%% check inputs
est = 'trimmean';
trimming = 20/100; % amount of trimming
plot_option = 1;
alphav = 5/100;
nboot = 600;
X = varargin{1};
Y = varargin{2};
    
if nargin>1 && nargin<7
    if nargin >= 3; est = varargin{3}; end
    EST{1} = 'mean'; EST{2} = 'median'; EST{3} = 'trimmean';
    if ~ischar(est)
        error('estimator argument must be a character')
    elseif isempty(strmatch(est,EST))
        error('unknown estimator name')
    end
    
    if nargin >= 4; plot_option = varargin{4}; end
    if nargin >= 5; alphav = varargin{5}; end
    if ~isnumeric(alphav)
        error('alphav must be a numeric')
    elseif alphav > 1
        alphav = alphav / 100;
    end
    
    if nargin == 6; nboot = varargin{6}; end
    if ~isnumeric(nboot)
        error('nboot must be a numeric')
    end
else
    error('wrong number of arguments')
end

% check data
if size(X) ~= size(Y)
    error('X and Y must be the same size')
end

[n,p] = size(X);
if n == 1 && p > 2
    X = X'; Y = Y'; % put row vectors as columns
end

% adjust p value if needed
alphav = alphav / p; % Bonferoni adjustment


%% compute the difference(s) and call rst_1ttest
data = X-Y;
if strcmp(est,'mean')
    [h,CI,p] = rst_1ttest(data,'mean',0,alphav,nboot);
elseif strcmp(est,'trimmean')
    [h,CI,p] = rst_1ttest(data,'trimmean',0,alphav,nboot);
else
    [h,CI,p] = rst_1ttest(data,'median',0,alphav,nboot);
end
 
%% plot
if plot_option == 1
    if strcmp(est,'mean')
        figure('Name', 'bootstrap paired t-test')
        ylabel('Mean difference','FontSize',14); hold on; grid on
        if p ==1
            title('Mean and 95% CI','FontSize',16); grid on
        else
            title(['Means and ' num2str(100-alphav*100) '% CI'],'FontSize',16); grid on
        end
        bar([1:size(data,2)],nanmean(data)); 
        errorbar([1:size(data,2)],nanmean(data),nanmean(data)-CI(1,:),nanmean(data)-CI(2,:),'r','LineWidth',2);
        
    else
        figure('Name', 'Percentile paired t-test')
        if strcmp(est,'median')
            ylabel('Median difference','FontSize',14); hold on; grid on
            if p ==1
                title('Median and 95% CI','FontSize',16); grid on
            else
                title(['Medians and ' num2str(100-alphav*100) '% CI'],'FontSize',16); grid on
            end
            bar([1:size(data,2)],nanmedian(data));
            errorbar([1:size(data,2)],nanmedian(data),nanmedian(data)-CI(1,:),nanmedian(data)-CI(2,:),'r','LineWidth',2);
             
        elseif strcmp(est,'trimmean')
            ylabel('TrimMean difference','FontSize',14); hold on; grid on
            if p ==1 
                title('TrimMean and 95% CI','FontSize',16); grid on
            else
                title(['TrimMeans and ' num2str(100-alphav*100) '% CI'],'FontSize',16); grid on
            end
            bar([1:size(data,2)],rst_trimmean(data));
            errorbar([1:size(data,2)],rst_trimmean(data),rst_trimmean(data)-CI(1,:),rst_trimmean(data)-CI(2,:),'r','LineWidth',2);
        end
    end
        
    if size(data,2) ==1
        xlabel('Pair 1','FontSize',14);        
    else
        for i=1:size(data,2)
            xname{i} = [num2str(i) ' '];
        end
        xlabel(sprintf('Pairs %s', cell2mat(xname)),'FontSize',14);
    end
    set(gca,'FontSize',12);
    v=axis; ymin = v(3) + (v(4)-v(3))/10;
    set(gcf,'Color','w')
end    