function [h,M,CI,p] = rst_2ttest(varargin)

% performs a robust two samplea test, ie a percentile bootstrap on the
% difference of the estimators
%
% FORMAT
% [h,M,CI] = rab_2ttest(X,Y)
% [h,M,CI,p] = rab_2ttest(X,Y,'estimator',plot_option,alpha,nboot)
%
% INPUTS
% - X and Y are a n-by-p matrix, if p>1 the test is computed for successive columns 
% and resampling is identical in each column
% - estimator = 'mean' or 'median' or 'trimmean'
% - plot_option indicates to plot the result (1 - default) or not (0)
% - alpha is the alpha level, default = 5%
% - nboot is the number of bootstrap, default = 600
%
% OUTPUTS
% [h,M,CI] = rab_1ttest(X,Y,'mean')
% h is the hypothesis that the mean is 0 if h = 0 then this is not
% significantly different from 0 ; note that if data is a matrix, 
% a Bonferonni correction is applied for the number of tests (columns)
% M is the mean difference
% CI are the alpha percent (Bonferoni corrected) confidence intervals used to get h
%
% [h,M,CI,p] = rab_1ttest(X,Y,'median');
% [h,M,CI,p] = rab_1ttest(X,Y,'trimmean');
% h is the hypothesis that the median or 20% trimmean is 0 if h = 0 then this is not
% significantly different from 0 ; note that if data is a matrix, 
% a Bonferonni correction is applied for the number of tests (columns)
% M is the central tendency estimator of the difference
% CI are the alpha percent (Bonferoni corrected) confidence intervals
% p are the p values
%
% Cyril Pernet - v1 - April 2012
% -------------------------------

%% check arguments
est = 'trimmean';
trimming = 20/100; % amount of trimming
alphav = 5/100;
nboot = 600;
plot_option = 1;
X = varargin{1};
Y = varargin{2};
 
[n,p] = size(X);
if n == 1 && p > 1 % if vectors are used, doesn't matter the direction
    X = X';
    Y = Y';
end

% other argument
if nargin>=2 && nargin<7
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

%% adjust p value if needed
alphav = alphav / size(X,2); % Bonferoni adjustment

%% create a resampling table

table1 = randi(size(X,1),size(X,1),nboot);
table2 = randi(size(Y,1),size(Y,1),nboot);
for c=1:size(X,2)
        v = X(:,c);
        data1{c} = v(table1); % applies the sample resampling for each column
        v = Y(:,c);
        data2{c} = v(table2); % applies the sample resampling for each column
end

%% do the percentile test

low = round((alphav*nboot)/2);
high = nboot - low;
n(1) = size(X,1);
n(2) = size(Y,1);

% get the difference just for output
if strcmp(est,'mean')
    M = nanmean(X,1) - nanmean(Y,1);
elseif strcmp(est,'median')
    M = nanmedian(X,1) - nanmedian(Y,1);
elseif strcmp(est,'trimmean')
    M = rab_trimmean(X,trimming)- rab_trimmean(Y,trimming);
end

for c=1:size(X,2)
    if strcmp(est,'mean')
        m{c} = sort(nanmean(data1{c},1) - nanmean(data2{c},1));
    elseif strcmp(est,'median')
        m{c} = sort(nanmedian(data1{c},1) - nanmedian(data2{c},1));
    elseif strcmp(est,'trimmean')
        m{c} = sort(rab_trimmean(data1{c},trimming)- rab_trimmean(data2{c},trimming));
    end
    
    CI(1,c) = m{c}(low);
    CI(2,c) = m{c}(high);
    
    pb = sum(m{c} > 0) / nboot;
    p(c) = 2*min(pb,1-pb);
end

h = p<alphav;

%% plot
    
if plot_option == 1
    figure('Name', 'Two-samples percentile bootstrapt')
    if strcmp(est,'mean')
        ylabel('Mean diffence','FontSize',14); hold on; grid on
        if size(X,2) ==1
            title('Mean and 95% CI','FontSize',16); grid on
        else
            title(['Means and ' num2str(100-alphav*100) '% CI'],'FontSize',16); grid on
        end
        
    elseif strcmp(est,'median')
        ylabel('Median difference','FontSize',14); hold on; grid on
        if size(X,2) ==1
            title('Median and 95% CI','FontSize',16); grid on
        else
            title(['Medians and ' num2str(100-alphav*100) '% CI'],'FontSize',16); grid on
        end
        
    elseif strcmp(est,'trimmean')
        ylabel('TrimMean difference','FontSize',14); hold on; grid on
        if size(X,2) ==1
            title('TrimMean and 95% CI','FontSize',16); grid on
        else
            title(['TrimMeans and ' num2str(100-alphav*100) '% CI'],'FontSize',16); grid on
        end
    end
    
    bar([1:size(X,2)],M);
    errorbar([1:size(X,2)],M,M-CI(1,:),M-CI(2,:),'r','LineWidth',2);
        
    if size(X,2) ==1
        xlabel('Condition 1','FontSize',14);        
    else
        for i=1:size(X,2)
            xname{i} = [num2str(i) ' '];
        end
        xlabel(sprintf('Conditions %s', cell2mat(xname)),'FontSize',14);
    end
    set(gca,'FontSize',12);
    v=axis; ymin = v(3) + (v(4)-v(3))/10;
    set(gcf,'Color','w')
end


