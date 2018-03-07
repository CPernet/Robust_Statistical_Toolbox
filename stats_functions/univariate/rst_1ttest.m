function [h,CI,p] = rst_1ttest(varargin)

% performs a robust one sample test, ie tests if an estimator is different
% from zeros.
%
% FORMAT
% [h,CI] = rst_1ttest(data)
% [h,CI,p] = rst_1ttest(data,'estimator',plot_option,alpha,nboot)
%
% INPUTS
% - data is either a vector or a matrix, in the latter case tests are
% performed column wise and the resampling is identical in each column
% - estimator = 'mean' --> performs a bootstrap-t on the mean
% - estimator = 'median' --> performs a percentile bootstrap on the median
%                note the median is the mid decile Harell-Davis estimator
% - estimator = 'trimmean' --> performs a percentile bootstrap on the 20% trimmean (default)
% - plot_option indicates to plot the result (1 - default) or not (0)
% - alpha is the alpha level, default = 5%
% - nboot is the number of bootstrap, default = 600
%
% OUTPUTS
% [h,CI] = rst_1ttest(data,'mean')
% h is the hypothesis that the mean is 0 if h = 0 then this is not
% significantly different from 0 ; note that if data is a matrix,
% a Bonferonni correction is applied for the number of tests (columns)
% CI are the alpha percent (Bonferoni corrected) confidence intervals used to get h
%
% [h,CI,p] = rst_1ttest(data,'median');
% [h,CI,p] = rst_1ttest(data,'trimmean');
% h is the hypothesis that the median or 20% trimmean is 0 if h = 0 then this is not
% significantly different from 0 ; note that if data is a matrix,
% a Bonferonni correction is applied for the number of tests (columns)
% p are the p values
% CI are the alpha percent (Bonferoni corrected) confidence intervals
%
% Cyril Pernet - v1 - March 2012
% ----------------------------------------------------------------------
% Copyright (C) RST Toolbox Team 2014


%% check arguments
est = 'trimmean';
trimming = 20/100; % amount of trimming
alphav = 5/100;
nboot = 600;
plot_option = 1;
data = varargin{1};

if nargin>=1 && nargin<6
    if nargin >= 2; est = varargin{2}; end
    EST{1} = 'mean'; EST{2} = 'median'; EST{3} = 'trimmean';
    if ~ischar(est)
        error('estimator argument must be a character')
    elseif isempty(strmatch(est,EST))
        error('unknown estimator name')
    end
    
    if nargin >= 3; plot_option = varargin{3}; end
    if nargin >= 4; alphav = varargin{4}; end
    if ~isnumeric(alphav)
        error('alphav must be a numeric')
    elseif alphav > 1
        alphav = alphav / 100;
    end
    
    if nargin == 5; nboot = varargin{5}; end
    if ~isnumeric(nboot)
        error('nboot must be a numeric')
    end
else
    error('wrong number of arguments')
end

% check data
[n,p] = size(data);
if n == 1 && p > 2
    data = data';
end
[n,p] = size(data);

% adjust p value if needed
alphav = alphav / p; % Bonferoni adjustment

%% create a resampling table

% resample the column(s)
table = randi(n,n,nboot);
if p == 1
    D = data(table);
else
    for c = 1:p
        v = data(:,c);
        D{c} = v(table); % applies the sample resampling for each column
    end
end


%% performs a bootstrap-t

if strcmp(est,'mean')
    
    % bootstrap parameters for symmetric CI
    boot_param = round((1-alphav)*nboot);
    
    % now down to the resampling
    if p==1
        try
            M = nanmean(D);
        catch ME
            error('NaN in the data, mean estimator cannot be computed without the stat toolbox or the NaN suite')
        end
        S = nanstd(D);
        if sum(abs(M-nanmean(data)))==0 && sum(S)==0
            T = zeros(1,length(S));
        else
            T = sort((sqrt(n).*abs(M-nanmean(data)))./S);
            T = T(isfinite(T)); boot_param = round((1-alphav)*length(T));
        end
        C = T(boot_param)*(nanstd(data)/sqrt(n));
        CI = [nanmean(data) - C; nanmean(data) + C];
        h  = (CI(1) > 0) + (CI(2) < 0); % lower CI above 0 or higher CI < 0
    else
        for c = 1:p
            try
                M = nanmean(D{c});
            catch ME
                error('NaN in the data, mean estimator cannot be computed without the stat toolbox or thew NaN suite')
            end
            S = nanstd(D{c});
            if sum(abs(M-nanmean(data(:,c))))==0 && sum(S)==0
                T = zeros(1,length(S));
            else
                T = sort((sqrt(n).*abs(M-nanmean(data(:,c))))./S);
                T = T(isfinite(T)); boot_param = round((1-alphav)*length(T));
            end
            C = T(boot_param)*(nanstd(data(:,c))/sqrt(n));
            CI(:,c) = [nanmean(data(:,c)) - C; nanmean(data(:,c)) + C];
            h(c)  = (CI(1,c) > 0) + (CI(2) < 0) ; % lower CI above 0
        end
    end
    
    % add an exeption
    for c=1:size(data,2)
        if sum(CI(:,c)) == 0
            h(c) = 0;
        end
    end
    
    % check arguments out
    if nargout == 3
        p = 'no p values for mean estimates';
    end
    
    %% performs a percentile bootstrap
else
    
    % bootstrap parameters
    low = round((alphav*nboot)/2);
    high = nboot - low;
    
    % now down to the resampling
    if strcmp(est,'median')
        try
            M = rst_hd(data,.5);
        catch ME
            error('NaN in the data, median estimator cannot be computed without the stat toolbox or the NaN suite')
        end
        if p==1
            Mb = sort(rst_hd(D,.5));
            if sum(isnan(Mb)) > 0
                Mb(find(isnan(Mb))) = [];
                low = round((alphav*length(Mb))/2);
                high = length(Mb) - low;
            end
            CI = [Mb(low+1) ; Mb(high)];
            pb = sum(Mb<0) / length(Mb);
        else
            for c = 1:p
                Mb = sort(rst_hd(D{c},.5));
                if sum(isnan(Mb)) > 0
                    Mb(find(isnan(Mb))) = [];
                    low = round((alphav*length(Mb))/2);
                    high = length(Mb) - low;
                end
                CI(:,c) = [Mb(low+1) ; Mb(high)];
                pb(c) = sum(Mb<0) / length(Mb);
            end
        end
    elseif strcmp(est,'trimmean')
        M = rst_trimmean(data,trimming); % default here as 20% trimmean
        % can only be changed by editing argument at the top
        if p==1
            Mb = sort(rst_trimmean(D,trimming));
            if sum(isnan(Mb)) > 0
                Mb(find(isnan(Mb))) = [];
                low = round((alphav*length(Mb))/2);
                high = length(Mb) - low;
            end
            CI = [Mb(low+1) ; Mb(high)];
            pb = sum(Mb<0) / length(Mb);
        else
            for c = 1:p
                Mb = sort(rst_trimmean(D{c},trimming));
                if sum(isnan(Mb)) > 0
                    Mb(find(isnan(Mb))) = [];
                    low = round((alphav*length(Mb))/2);
                    high = length(Mb) - low;
                end
                CI(:,c) = [Mb(low+1) ; Mb(high)];
                pb(c) = sum(Mb<0) / length(Mb);
            end
        end
    end
    
    % check which direction is the effect
    % add an execption
    for c=1:size(data,2)
        if sum(CI(:,c)) == 0
            p(c) = 1; h(c) = 0;
        else
            p(c) = min(pb(c), 1-pb(c));
            h(c) = p(c)<= (alphav/2);
            if p(c) == 0
                p(c) = 1/nboot;
            end
        end
    end
end

h = logical(h);

if plot_option == 1
    L = size(data,2);
    color_scheme = cubehelixmap('semi_continuous',L+10);
    
    if strcmp(est,'mean')
        figure('Name', 'One-sample bootstrap t-test')
        ylabel('Mean value','FontSize',14); hold on; grid on
        if p ==1
            title('Mean and 95% CI','FontSize',16); grid on
        else
            title(['Means and ' num2str(100-alphav*100) '% CI'],'FontSize',16); grid on
        end
        
        for i=1:L
            tmp = repmat(data(:,i),1,2); tmp(1:2:end,1) = NaN; tmp(2:2:end,2) = NaN;
            scatter(repmat(i-0.05,n,1),tmp(:,1),50,'k'); hold on
            scatter(repmat(i+0.05,n,1),tmp(:,2),50,'k');
            rectangle('Position',[i-0.2,CI(1,i),0.4,CI(2,i)-CI(1,i)],'Curvature',[0.4 0.4],'LineWidth',2,...
                'FaceColor',color_scheme(i+8,:),'EdgeColor',[0.35 0.35 0.35]);
            plot([i-0.2 i+0.2],[nanmean(data(:,i)) nanmean(data(:,i))],'LineWidth',3,'Color',[0.35 0.35 0.35]);
        end
        
        
    else
        figure('Name', 'One-sample percentile t-test')
        if strcmp(est,'median')
            ylabel('Median value','FontSize',14); hold on; grid on
            if p ==1
                title('Median and 95% CI','FontSize',16); grid on
            else
                title(['Medians and ' num2str(100-alphav*100) '% CI'],'FontSize',16); grid on
            end
            
            for i=1:L
                tmp = repmat(data(:,i),1,2); tmp(1:2:end,1) = NaN; tmp(2:2:end,2) = NaN;
                scatter(repmat(i-0.05,n,1),tmp(:,1),50,'k'); hold on
                scatter(repmat(i+0.05,n,1),tmp(:,2),50,'k');
                rectangle('Position',[i-0.2,CI(1,i),0.4,CI(2,i)-CI(1,i)],'Curvature',[0.4 0.4],'LineWidth',2,...
                    'FaceColor',color_scheme(i+8,:),'EdgeColor',[0.35 0.35 0.35]);
                plot([i-0.2 i+0.2],[nanmedian(data(:,i)) nanmedian(data(:,i))],'LineWidth',3,'Color',[0.35 0.35 0.35]);
            end
            
        elseif strcmp(est,'trimmean')
            ylabel('TrimMean value','FontSize',14); hold on; grid on
            if p ==1
                title('TrimMean and 95% CI','FontSize',16); grid on
            else
                title(['TrimMeans and ' num2str(100-alphav*100) '% CI'],'FontSize',16); grid on
            end
            
            for i=1:L
                tmp = repmat(data(:,i),1,2); tmp(1:2:end,1) = NaN; tmp(2:2:end,2) = NaN;
                scatter(repmat(i-0.05,n,1),tmp(:,1),50,'k'); hold on
                scatter(repmat(i+0.05,n,1),tmp(:,2),50,'k');
                rectangle('Position',[i-0.2,CI(1,i),0.4,CI(2,i)-CI(1,i)],'Curvature',[0.4 0.4],'LineWidth',2,...
                    'FaceColor',color_scheme(i+8,:),'EdgeColor',[0.35 0.35 0.35]);
                plot([i-0.2 i+0.2],[rst_trimmean(data(:,i)) rst_trimmean(data(:,i))],'LineWidth',3,'Color',[0.35 0.35 0.35]);
            end
            
            
        end
    end
    
    if size(data,2) ==1
        xlabel('Condition 1','FontSize',14);
    else
        for i=1:size(data,2)
            xname{i} = [num2str(i) ' '];
        end
        xlabel(sprintf('Conditions %s', cell2mat(xname)),'FontSize',14);
    end
    set(gca,'FontSize',12);
    v=axis; ymin = v(3) + (v(4)-v(3))/10;
    set(gcf,'Color','w')
end


