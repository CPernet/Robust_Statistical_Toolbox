function [est,HDI,OUT,K] = rst_data_plot(Data,varargin)

% plots the data split by groups showing each data points, marking S outliers,
% with a summary statistics estimator and it' 95&% Bayes bootstrap Highest Density Interval
% along with the the kernel density estimate of the data distribution (and
% interquartile range)
%
% FORMAT: [estimators,HDI,outliers,Kernel_densities] = rst_data_plot(Data,options)
%
% INPUT: Data is a matrix, data are plotted colmun-wise (= different groups to plot)
%        options are
%                'datascatter' can be 'on' (default) or 'off' to have data scatter plots for group
%                'within' a value for the distance/width between points in a group a.k.a. jitter
%                'point_size' a value for the size of data points
%                'between' a value for the distance between groups
%                'bubble' can be 'on' or 'off' (default) to change size based on the number of similar data points 
%                'outliers' can be 'on' (default) or 'off' to plot differently S-outliers 
%                'bars' can be 'on' or 'off' (default) to plot as bar graph (datascatter is then on) 
%                'estimator' can be 'median' (default), 'mean', 'trimmed mean'
%                'kernel' can be 'on' (default) or 'off' to add distribution
%                'kernel_plot' can be 'full' (draw doubled distributions on top of scatter plots),
%                                     'half' (scatter plots and distributions side by side)
%                                     'side' (scatter plots and distributions on the same x axis on the side)
%                'newfig' 'yes' or 'no' to plot in a new figure or the
%                          current one (useful for making subplots yourself)
%
% OUTPUT: estimators is a vector the summary statistics
%         HDI the 95% High Density Interval (Bayesian bootstrap)
%         Outliers is a binary matrix indicating S-outliers
%         Kernel_densities the kernel density estimated (returns bins and frequencies)
%
% see also rst_outlier, rst_RASH, rst_trimmean, rst_hd, rst_colour_map
%
% Examples:
% ---------
% show mean and 95% HDI as bars, all data, do not mark outliers, 
% distributions on the side with interquartile range
% Data = [randg(1,100,1) randn(100,1) randn(100,1) 1-randg(1,100,1)];
% Data(50:55,1) = NaN; Data(90:100,2) = NaN;
% [est,HDI,outliers,KDE]=rst_data_plot(Data,'outliers','off','estimator','mean');
%
% show trimmed mean and 95% HDI as a box, all data with outliers marked, 
% distribution on top (no interquartile range)
% Data = [randg(1,100,1) randn(100,1) randn(100,1) 1-randg(1,100,1)];
% Data(50:55,1) = NaN; Data(90:100,2) = NaN;
% [est,HDI,outliers,KDE]=rst_data_plot(Data,'estimator','trimmed mean','kernel_plot','full');
%
% show median and 95% HDI as bar plot, outliers only
% distributions on the far end with interquartile range
% Data = [randn(100,1) randg(1,100,1)];
% Data(50:55,1) = NaN; Data(90:100,2) = NaN;
% [est,HDI,outliers,KDE]=rst_data_plot(Data,'estimator','median','datascatter','off','point_size',50,'kernel_plot','side','bars','on');
%
% show median and 95% HDI with bubbles on to show repeted observations
% Data = randi(30,1,200)';
% [est,HDI,outliers,KDE]=rst_data_plot(Data,'estimator','median','bubble','on');
%
% Cyril Pernet
% Eric Nicholas (fixed threshold for similar data points)
% -------------------------------------------------------------------------
% Copyright (C) RST Toolbox Team 2021

%% Defaults

% hard coded
Nb = 1000;               % number of bootstrap samples
prob_coverage  = 95/100; % prob coverage of the HDI
outlier_method = 3;      % 1 MAD-median rule, 2 small sample adjusted, 3 S-outliers
dist_method    = 'ASH';  % Average Shifted Histogram, could also be RASH
trimming       = 20;     % 20% trimming each side of the distribution
decile         = .5;     % Median estimated using the 5th decile of Harell Davis

% soft coded, see options
estimator             = 'median';
datascatter           = 'on';
bars                  = 'off';
bubble                = 'off';
outlier_plot          = 'on';
between_gp_dispersion = 0.25;
within_gp_dispersion  = 0.025;
point_size            = 25;
kernel                = 'on';
kernel_plot           = 'half'; 
newfig                = 'yes';

%% check inputs
if ~exist('Data','var')
    [name,place,sts]=uigetfile({'*.csv;*.tsv;*.txt'},'Select data file');
    if sts == 0
        return
    else
        if strcmp(name(end-2:end),'csv')
            Data = csvread([place name]); % expect a comma between variables
        elseif strcmp(name(end-2:end),'tsv')
            Data = dlmread([place name]); % expect a tab between variables
        elseif strcmp(name(end-2:end),'txt')
            Data = load([place name]);    % assumes just a space between variables
        end
    end
end

if nargin ==1
    estimator = questdlg('Which summary statistics to plot?','Stat question','Mean',' 20% Trimmed mean','Median','Median');
    if strcmp(estimator,' 20% Trimmed mean')
        estimator = 'Trimmed mean';
    end
else
    for n=1:(nargin-1)
        if strcmpi(varargin(n),'between')
            between_gp_dispersion = cell2mat(varargin(n+1));
        elseif strcmpi(varargin(n),'within')
            within_gp_dispersion = cell2mat(varargin(n+1));
        elseif any(strcmpi(varargin(n),{'outlier','outliers'}))
            outlier_plot = cell2mat(varargin(n+1));
        elseif strcmpi(varargin(n),'point_size')
            point_size = cell2mat(varargin(n+1));
        elseif any(strcmpi(varargin(n),{'bar','bars'}))
            bars = cell2mat(varargin(n+1));
        elseif strcmpi(varargin(n),'estimator')
            if ~strcmpi(varargin(n+1),'median') && ...
                    ~strcmpi(varargin(n+1),'mean') && ...
                    ~strcmpi(varargin(n+1),'trimmed mean')
                error(['estimator ' cell2mat(varargin(n+1)) ' is not recognized'])
            else
                estimator = cell2mat(varargin(n+1));
            end
        elseif strcmpi(varargin(n),'kernel')
            kernel = cell2mat(varargin(n+1)); K = [];
        elseif any(strcmpi(varargin(n),{'kernel_plot','kernel plot'}))
            kernel_plot= cell2mat(varargin(n+1));
        elseif any(strcmpi(varargin(n),{'datascatter','data scatter'}))
            datascatter = cell2mat(varargin(n+1));
        elseif strcmpi(varargin(n),'bubble')
            bubble = cell2mat(varargin(n+1));
        elseif strcmpi(varargin(n),'newfig')
            newfig = cell2mat(varargin(n+1));
        end
    end
end

% with doubled kernel, some defaults are too small
if strcmpi(kernel_plot,'full') && within_gp_dispersion == 0.025
     within_gp_dispersion = 0.05;
end

if strcmpi(kernel_plot,'full') && point_size == 25
    point_size = 40;
end

% that bahaviour is not allowed 
if strcmpi(bars,'on') && strcmpi(datascatter,'off')
    warning('bar plots hide data, scatter turned on')
    datascatter = 'on';
end

%% how many groups, set gp_index (the x-axis ticks)
[N,grouping] = size(Data);
if N == 1 && grouping > 1
    Data = Data';
    [~,grouping] = size(Data);
end

if strcmpi(kernel,'on')
    if strcmpi(kernel_plot,'side')
        gp_index = linspace(1,1+between_gp_dispersion*grouping,grouping);
        gp_index = [gp_index max(gp_index)+0.5];
    else
        gp_index = 1:(1+between_gp_dispersion):(grouping*(1+between_gp_dispersion));
    end
else
    gp_index = linspace(1,1+between_gp_dispersion*grouping,grouping);
end

%% Summary stat
if strcmpi(estimator,'Mean')
    est = nanmean(Data,1);
elseif strcmpi(estimator,'Trimmed mean')
    est = rst_trimmean(Data,trimming);
elseif strcmpi(estimator,'Median')
    est = rst_hd(Data,decile); % Median estimated using the 5th decile of Harell Davis
end
HDI = zeros(2,grouping); % placeholder for CI

%% start
if strcmpi(newfig,'yes')
    figure; 
end
hold on

%% select color scheme
color_scheme = flipud(rst_colour_maps(grouping));

%% Scatter plot of the data with automatic spread
% scatter plot parameters

if strcmp(kernel,'on')
    disp('Computing distribution kernels and plotting ...')
else
    disp('start plotting ...')
end

OUT = Data;
cst = max(abs(diff(Data(:)))) * 0.1;
for u=1:grouping
    if strcmp(bubble,'on')
        tmp  = sort(Data(~isnan(Data(:,u)),u)); 
    else
        tmp  = Data(~isnan(Data(:,u)),u); 
    end
    outliers = rst_outlier(tmp,outlier_method); % find outliers
    OUT(~isnan(Data(:,u)),u) = outliers; 
    OUT = logical(OUT);

    if strcmpi(datascatter,'on') || strcmpi(outlier_plot,'on') % plot individual data points
        
        % look for overlapping data points in the display
        if length(unique(tmp)) == length(tmp)
            change = diff(tmp) < round(range(tmp)/point_size);
        else
            change = ~diff(tmp);
        end
        
        if isempty(find(change,1))
            bubble = 'off'; 
        end
        
        if strcmp(bubble,'off')
            Y = [tmp tmp];                                     % make it two columns to display 
            S = repmat(point_size,[length(tmp),1]);            % constant point size
            Y(1:2:length(tmp),1) = NaN;                        % plot every other point 
            Y(2:2:length(tmp),2) = NaN;
            outliers = single([outliers outliers]);            % same for outliers
            outliers(1:2:length(tmp),1) = NaN; 
            outliers(2:2:length(tmp),2) = NaN;
        else
            tmp(logical([0;change])) = [];                     % remove overlapping data points (since they are indicated by size)
            Y = [tmp tmp];                                     % make it two columns to display 
            Y(1:2:length(tmp),1) = NaN;                        % plot every other point
            Y(2:2:length(tmp),2) = NaN;
            outliers(logical([0;change])) = [];                % same for outliers
            outliers = single([outliers outliers]);
            outliers(1:2:length(outliers),1) = NaN; 
            outliers(2:2:length(outliers),2) = NaN;
            for c=length(tmp):-1:1
                original_data = Data(~isnan(Data(:,u)),u);
                S(c) = sum(original_data == tmp(c))*point_size; % S is the size of the data points
            end
            S = S-min(S)+point_size;
        end
        
        % plot
        % ----
        X = repmat([-within_gp_dispersion/2 within_gp_dispersion/2],[length(tmp),1]) + gp_index(u);
        if strcmpi(kernel,'on') && strcmpi(kernel_plot,'half')
            X = X - 0.25; % shift left
        end
        
        if strcmpi(bars,'on')
            h = bar(mean(X(1,:)),est(u),within_gp_dispersion*3);
            set(h,'BaseValue',min(Data(:))-cst,'FaceAlpha',0.5,'FaceColor',color_scheme(u,:),'EdgeColor',color_scheme(u,:)); hold on
        end
        
        for p=1:size(Y,2)
            if strcmpi(datascatter,'on')
                if strcmp(bubble,'off') || isempty(find(change))
                    if strcmpi(outlier_plot,'on')
                        if strcmpi(bars,'on')
                            scatter(X(outliers(:,p)==0,p),Y(outliers(:,p)==0,p),S(find(outliers(:,p)==0)),[0.5 0.5 0.5]); hold on;
                        else
                            scatter(X(outliers(:,p)==0,p),Y(outliers(:,p)==0,p),S(find(outliers(:,p)==0)),color_scheme(u,:)); hold on;
                        end
                    else
                        if strcmpi(bars,'on')
                            scatter(X(:,p),Y(:,p),S,[0.5 0.5 0.5]); hold on;
                        else
                            scatter(X(:,p),Y(:,p),S,color_scheme(u,:)); hold on;
                        end
                    end
                else % bubble is on, fill the markers
                    if strcmpi(outlier_plot,'on')
                        if strcmpi(bars,'on')
                            scatter(X(outliers(:,p)==0,p),Y(outliers(:,p)==0,p),S(find(outliers(:,p)==0)),[0.5 0.5 0.5],'filled'); hold on;
                        else
                            scatter(X(outliers(:,p)==0,p),Y(outliers(:,p)==0,p),S(find(outliers(:,p)==0)),color_scheme(u,:),'filled'); hold on;
                        end
                    else
                        if strcmpi(bars,'on')
                            scatter(X(:,p),Y(:,p),S,[0.5 0.5 0.5],'filled'); hold on;
                        else
                            scatter(X(:,p),Y(:,p),S,color_scheme(u,:),'filled'); hold on;
                        end
                    end
                end
            end
            
            if strcmpi(outlier_plot,'on')
                if strcmpi(kernel_plot,'full') % make sure it's visible
                    if strcmpi(bars,'on')
                        scatter(X(outliers(:,p)==1,p),Y(outliers(:,p)==1,p),S(find(outliers(:,p)==1)).*2,[0.5 0.5 0.5],'+');
                    else
                        scatter(X(outliers(:,p)==1,p),Y(outliers(:,p)==1,p),S(find(outliers(:,p)==1)).*2,color_scheme(u,:),'+');
                    end
                else
                    if strcmpi(bars,'on')
                        scatter(X(outliers(:,p)==1,p),Y(outliers(:,p)==1,p),(S(find(outliers(:,p)==1))-(point_size/2)),[0.5 0.5 0.5],'+');
                    else
                        scatter(X(outliers(:,p)==1,p),Y(outliers(:,p)==1,p),(S(find(outliers(:,p)==1))-(point_size/2)),color_scheme(u,:),'+');
                    end
                end
            end
        end
    end % end of scatter plot
     
    %% Add the density estimate
    if strcmp(kernel,'on') && sum(tmp) ~=0
        % get the kernel
        [bc,K]=rst_RASH(tmp,100,dist_method);
        % remove 0s
        bc(K==0)=[]; K(K==0)= [];
        % create symmetric values
        K = (K - min(K)) ./ max(K); % normalize to [0 1] interval
        high=(K/2); low=(-high);
        
        % plot contours
        if strcmpi(kernel_plot,'side')
            y1 = plot(high+gp_index(end),bc);
        else
            y1 = plot(high+gp_index(u),bc);
        end
        set(y1,'Color',color_scheme(u,:)); hold on
        if isnumeric(y1)
            y1 = get(y1); % for older matlab versions
        end
        
        if strcmpi(kernel_plot,'full')
            y2 = plot(low+gp_index(u),bc); 
            set(y2,'Color',color_scheme(u,:));
            if isnumeric(y2)
                y2 = get(y2); 
            end
            xpoints = [y2.XData',y1.XData']; % fill
            filled  = [y2.YData',y1.YData'];
        else
            xpoints = [repmat(min(y1.XData),1,length(y1.XData))',y1.XData']; 
            filled  = [y1.YData',y1.YData'];            
        end
        
        % check that we have continuous data, otherwise 0 pad
        if diff(xpoints(1,:)) > 0.001*(range(xpoints(:,1)))
            xtmp          = NaN(size(xpoints,1)+1,size(xpoints,2));
            if strcmpi(kernel_plot,'side')
                xtmp(1,:) = [gp_index(end) gp_index(end)];
            else
                xtmp(1,:) = [gp_index(u) gp_index(u)];
            end
            xtmp(2:end,:) = xpoints;
            xpoints       = xtmp; clear xtmp
            ytmp          = NaN(size(filled,1)+1,size(filled,2));
            ytmp(1,:)     = filled(1,:);
            ytmp(2:end,:) = filled;
            filled        = ytmp; clear ytmp
        end
        
        if diff(xpoints(end,:)) > 0.001*(range(xpoints(:,1)))
            if strcmpi(kernel_plot,'side')
                xpoints(end+1,:) = [gp_index(end) gp_index(end)];
                filled(end+1,:)  = filled(end,:);
            else
                xpoints(end+1,:) = [gp_index(u) gp_index(u)];
                filled(end+1,:)  = filled(end,:);
            end
        end
        hold on; 
        fillhandle = fill(xpoints,filled,color_scheme(u,:));
        set(fillhandle,'LineWidth',2,'EdgeColor',color_scheme(u,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
        
        %% add IQR - using again Harell-Davis Q
        if any(strcmpi(kernel_plot,{'half','side'}))
            % lower quartile
            ql           = rst_hd(tmp,0.25);
            [~,position] = min(abs(filled(:,1) - ql));
            plot(xpoints(position,:),filled(position,:),'Color',color_scheme(u,:));
            
            % upper quartile
            qu = rst_hd(tmp,0.75);
            [~,position] = min(abs(filled(:,1) - qu));
            plot(xpoints(position,:),filled(position,:),'Color',color_scheme(u,:));
            
            % median
            if  strcmpi(estimator,'median')
                qm = est(u);
            else
                qm= rst_hd(tmp,0.5);
            end
            [~,position] = min(abs(filled(:,1) - qm));
            plot(xpoints(position,:),filled(position,:),'Color',color_scheme(u,:));
            clear xpoints filled
        end
    end
    
    %% Bayes bootstrap
    if Nb ~= 0
        % sample with replacement from Dirichlet
        % sampling = number of observations
        tmp = sort(Data(~isnan(Data(:,u)),u));
        n   =size(tmp,1); bb = zeros(Nb,1);
        fprintf('Computing Bayesian Bootstrap for HDI in group %g/%g\n',u,grouping)
        for boot=1:Nb % bootstrap loop
            theta    = exprnd(1,[n,1]);
            weigths  = theta ./ repmat(sum(theta),n,1);
            resample = (datasample(tmp,n,'Replace',true,'Weights',weigths));

            % compute the estimator
            if strcmpi(estimator,'Mean')
                bb(boot) = nanmean(resample,1);
            elseif strcmpi(estimator,'Trimmed Mean')
                bb(boot) = rst_trimmean(resample,trimming);
            elseif strcmpi(estimator,'Median')
                bb(boot) = rst_hd(resample,decile);
            end
        end
        sorted_data   = sort(bb); % sort bootstrap estimates
        upper_centile = floor(prob_coverage*size(sorted_data,1)); % upper bound
        nCIs          = size(sorted_data,1) - upper_centile;
        ci            = 1:nCIs;
        ciWidth       = sorted_data(ci+upper_centile) - sorted_data(ci); % all centile distances
        [~,index]     = find(ciWidth == min(ciWidth)); % densest centile
        if length(index) > 1 % many similar values
            index = index(1);
        end
        HDI(1,u)      = sorted_data(index);
        HDI(2,u)      = sorted_data(index+upper_centile);

        % plot this with error bars
        X = linspace(X(1,1)-within_gp_dispersion, X(1,2)+within_gp_dispersion, 7);
        if strcmpi(bars,'off')
            if strcmpi(bubble,'on')
                plot(X,repmat(est(u),[1,length(X)]),'LineWidth',4,'Color',[0.5 0.5 0.5]);
            else
                plot(X,repmat(est(u),[1,length(X)]),'LineWidth',4,'Color',color_scheme(u,:));
            end
        end

        if strcmpi(kernel_plot,'full')
            rectangle('Position',[X(1),HDI(1,u),X(7)-X(1),HDI(2,u)-HDI(1,u)],'Curvature',[0.4 0.4],'LineWidth',3,'EdgeColor',color_scheme(u,:))
        else
            if strcmpi(bars,'on')
                errorbar(X(4),est(u),est(u)-HDI(1,u),HDI(2,u)-est(u),'LineWidth',3,'Color',[0.5 0.5 0.5]);
            else
                errorbar(X(4),est(u),est(u)-HDI(1,u),HDI(2,u)-est(u),'CapSize',length(X)*2,'LineWidth',3,'Color',[0.5 0.5 0.5]);
            end
        end
    end
end

%% finish off
plot_min = min(Data(:))-cst;
if min(HDI(1,:)) < plot_min
    plot_min = min(HDI(1,:));
end
plot_max = max(Data(:))+cst;
if max(HDI(2,:)) > plot_max
    plot_max = max(HDI(2,:));
end
axis([0.3 gp_index(end)+0.7 plot_min plot_max])

if size(Data,1) == 1
    if strcmp(datascatter,'on') && strcmp(kernel,'off')
        title(sprintf('Data scatter with %s \n and 95%% High Density Interval',estimator),'FontSize',16);
    end
    if strcmp(kernel,'on')
        title(sprintf('Data distribution with %s \n and 95%% High Density Interval',estimator),'FontSize',16);
    end
else
    if strcmp(datascatter,'on') && strcmp(kernel,'off')
        title(sprintf('Data scatters with %s \n and 95%% High Density Interval',estimator),'FontSize',16);
    end
    if strcmp(kernel,'on')
        title(sprintf('Data distributions with %ss \n and 95%% High Density Intervals',estimator),'FontSize',16);
    end
end
grid on; box on; drawnow

% add an output if not specified during the call
if nargout == 0
    S = questdlg('Save computed data?','Save option','Yes','No','No');
    if strcmp(S,'Yes')
        if exist(place,'var')
            place = pwd;
        end
        cd(place); tmp = [est; HDI];
        if strcmpi(estimator,'Mean')
            save('Mean_and_HDI.txt','tmp','-ascii');
        elseif strcmpi(estimator,'Trimmed Mean')
            save('Trimmed-Mean_and_HDI.txt','tmp','-ascii');
        elseif strcmpi(estimator,'Median')
            save('Median_and_HDI.txt','tmp','-ascii');
        end
        save('S-outliers.txt','Outliers','-ascii');
        fprintf('Data saved in %s\n',place)
    end
end

% if exist('plotly','file') == 2
%     output = questdlg('Do you want to output this graph with Plotly?', 'Plotly option');
%     if strcmp(output,'Yes')
%         save_dir = uigetdir('select directory to save html files','save in');
%         if ~isempty(save_dir)
%             cd(save_dir)
%             try
%                 fig2plotly(gcf,'strip',true,'offline', true);
%             catch
%                 fig2plotly(gcf,'strip',true); % in case local lib not available
%             end
%         else
%             return
%         end
%     else
%         return
%     end
% end

