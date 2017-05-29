function [h,M,CI,p] = rst_2ttest(gp1,gp2,varargin)

% performs a robust two samples test, ie a percentile bootstrap on the
% difference of the estimators
%
% FORMAT
% [h,M,CI]   = rst_2ttest(X,Y)
% [h,M,CI,p] = rst_2ttest(X,Y,'alphav',0.05,'estimator','median','newfig','yes')
%
% INPUTS
% - X and Y are a n-by-p matrix, if p>1 the test is computed for successive columns
%   and resampling is identical in each column
% - alphav is the alpha level, default = 5%
% - 'estimator' can be 'mean', 'median' or 'trimmed mean'
% - newfig    = 'yes' or 'no' allows to plot within an existing one
%
% OUTPUTS
% h  is the hypothesis that the mean is 0 if h = 0 then this is not
%    significantly different from 0 ; note that if data is a matrix,
%    a Bonferonni correction is applied for the number of tests (columns)
% M  is the estimator difference
% CI is the alpha percent (Bonferoni corrected) confidence intervals used to get h
%
% Cyril Pernet - v2 - March 2017
% ----------------------------------------------------------------------
% Copyright (C) RST Toolbox Team 2017

%% check arguments

% _hard coded_
nboot = 600; % how many bootstraps are used

% _soft coded_
alphav   = 5/100; % non adjusted alpha value
estimator ='trimmed mean'; % by default use a robust estimator
newfig    = 'yes'; % start a new figure
plothist  = 'no'; % do not create additional histograms of differences

% _check inputs_
if any(size(gp1,2) ~= size(gp2,2))
    error('input data must be of the same size')
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

[n,p] = size(gp1);
if n == 1 && p > 1 % if vectors are used, doesn't matter the direction
    gp1 = gp1';
    gp2 = gp2';
end

%% adjust p value if needed
alphav = alphav / size(gp1,2); % Bonferoni adjustment
low = round(nboot.*orig_alpha./2); % CI bounds
high = nboot - low;

%% create a resampling table

table1 = randi(size(gp1,1),size(gp1,1),nboot);
table2 = randi(size(gp2,1),size(gp2,1),nboot);
for c=1:size(gp1,2)
    v = gp1(:,c);
    data1{c} = v(table1); % applies the sample resampling for each column
    v = gp2(:,c);
    data2{c} = v(table2); % applies the sample resampling for each column
end

%% do the percentile test

n(1) = size(gp1,1);
n(2) = size(gp2,1);

% get the difference just for output
if strcmpi(estimator,'mean')
    M = nanmean(gp1,1) - nanmean(gp2,1);
elseif strcmpi(estimator,'median')
    M = rst_hd(gp1,0.5) - rst_hd(gp2,0.5);
elseif strcmpi(estimator,'trimmed mean')
    M = rst_trimmean(gp1)- rst_trimmean(gp2);
end

for c=1:size(gp1,2)
    if strcmpi(estimator,'mean')
        m{c} = sort(nanmean(data1{c},1) - nanmean(data2{c},1));
    elseif strcmpi(estimator,'median')
        m{c} = sort(rst_hd(data1{c}) - rst_hd(data2{c}));
    elseif strcmpi(estimator,'trimmed mean')
        m{c} = sort(rst_trimmean(data1{c})- rst_trimmean(data2{c}));
    end
    
    CI(1,c) = m{c}(low);
    CI(2,c) = m{c}(high);
    
    pb = sum(m{c} > 0) / nboot;
    p(c) = 2*min(pb,1-pb);
    if p(c) == 0
       p(c) = 1/nboot; 
    end
end

h = p<alphav;

%% plot

if strcmpi(newfig,'yes')
    figure('Name', 'Two-samples percentile bootstrapt')
    set(gcf,'Color','w');
end

L = size(gp1,2);
color_scheme = cubehelixmap('semi_continuous',L+10);
if L ==1
    title([estimator ' and 95% CI'],'FontSize',16); grid on
else
    title([estimator 's and ' num2str(100-alphav*100) '% CI'],'FontSize',16); grid on
end
ylabel([estimator ' difference'],'FontSize',14); hold on; grid on

for i=1:L
    rectangle('Position',[i-0.2,CI(1,i),0.4,CI(2,i)-CI(1,i)],'Curvature',[0.4 0.4],'LineWidth',2,...
        'FaceColor',color_scheme(i+8,:),'EdgeColor',[0.35 0.35 0.35]); hold on;
    plot([i-0.2 i+0.2],[M(i) M(i)],'LineWidth',3,'Color',[0.35 0.35 0.35]);
    tmp = repmat(m{i}',1,2); tmp(1:2:end,1) = NaN; tmp(2:2:end,2) = NaN;
    scatter(repmat(i-0.025,nboot,1),tmp(:,1),50,'k');
    scatter(repmat(i+0.025,nboot,1),tmp(:,2),50,'k');
end

if L ==1
    xlabel('Condition 1','FontSize',14);
else
    for i=1:size(gp1,2)
        xname{i} = [num2str(i) ' '];
    end
    xlabel(sprintf('Differences %s', cell2mat(xname)),'FontSize',14);
end
v=axis; ymin = v(3) + (v(4)-v(3))/10;

if strcmpi(plothist,'yes')
    figure('Name', 'Histograms of differences')
    set(gcf,'Color','w'); L = size(gp1,2);
    nrow = round(L/5);
    for i=1:L
        subplot(nrow,5,i); rst_density_hist(m{i});
    end
end

