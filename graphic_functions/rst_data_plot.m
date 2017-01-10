function [est,HDI]=rst_data_plot_fix(Data,varargin)

% plots the data split by groups showing each data points with the
% distribution and a summary statistics estimator with 95&% Bayes boot HDI
%
% FORMAT: [est,HDI]= rst_data_plot(Data,options)
%         [est,HDI]= rst_data_plot(Data,'between',0.25,'within',0.025,'pointsize',50,'estimator','median') 
%
% INPUT: Data is a matrix, data are plotted colmun-wise
%        options are 
%                'between' for the distance between distributions
%                'within' for the distance/width between points in a group
%                'point_size' the size of data points
%                'estimator' can be 'median', 'mean', 'trimmed mean' 
% OUTPUT: est is the summary statistics
%         HDI the 95% High Density Interval (Bayes bootstrap)
%         Ouliers is a binary matrix indicating S-outliers
%         KDE is the hernel density estimated (returns bins and frequencies)
%
% example: Data = [randn(100,1) randg(1,100,1) randn(100,1) 1-randg(1,100,1)];
%          Data(50:55,1) = NaN; Data(90:100,2) = NaN;
%          [est,HDI]=rst_data_plot(Data,'between',0.25,'within',0.025,'pointsize',50,'estimator','median')
%
% see also cubehelixmap, rst_outlier, rst_RASH, rst_trimmean, rst_hd
%
% Cyril Pernet - The University of Edinburgh
% Eric Nicholas - The University of Rochester
% -------------------------------------------------------------------------
% Copyright (C) RST Toolbox Team 2015

%% Defaults

% hard coded
Nb = 1000;              % number of bootstrap samples 
prob_coverage = 95/100; % prob coverage of the HDI
outlier_method = 3;     % 1 MAD-median rule, 2 small sample adjusted, 3 S-outliers
dist_method = 'RASH';   % Random Average Shifted Histogram, could also be ASH
trimming = 20;          % 20% trimming each side of the distribution
decile = .5;            % Median estimated using the 5th decile of Harell Davis 

% soft coded, see options
between_gp_dispersion = 0.25;
within_gp_dispersion = 0.025;
point_size = 50;

% check inputs
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
            Data = load([place name]); % assumes just a space between variables
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
        elseif strcmpi(varargin(n),'point_size')
            point_size = cell2mat(varargin(n+1));
        elseif strcmpi(varargin(n),'estimator')
            if ~strcmpi(varargin(n+1),'median') && ...
                    ~strcmpi(varargin(n+1),'mean') && ...
                    ~strcmpi(varargin(n+1),'trimmed mean')
                error(['estimator ' cell2mat(varargin(n+1)) ' is not recognized'])
            else
                estimator = cell2mat(varargin(n+1));
            end
        end
    end
end

%% how many groups
grouping = size(Data,2);
gp_index = 1:(1+between_gp_dispersion):(grouping*(1+between_gp_dispersion));

%% Summary stat 
if strcmpi(estimator,'Mean')
    est = nanmean(Data,1);
elseif strcmpi(estimator,'Trimmed mean')
    est = rst_trimmean(Data,trimming); 
elseif strcmpi(estimator,'Median')
    est = rst_hd(Data,decile); % Median estimated using the 5th decile of Harell Davis
end

%% Bayes bootstrap

% sample with replcaement from Dirichlet
% sampling = number of observations, e.g. participants
n=size(Data,1); 
bb = zeros(Nb,grouping);
for boot=1:Nb % bootstrap loop
    theta = exprnd(1,[n,1]);
    weigths = theta ./ repmat(sum(theta),n,1);
    resample = (datasample(Data,n,'Replace',true,'Weights',weigths));
    
    % compute the estimator
    if strcmpi(estimator,'Mean')
        bb(boot,:) = nanmean(resample,1);
    elseif strcmpi(estimator,'Trimmed Mean')
        bb(boot,:) = rst_trimmean(resample,20); 
    elseif strcmpi(estimator,'Median')
        bb(boot,:) = rst_hd(resample,.5); 
    end
end

sorted_data = sort(bb,1); % sort bootstrap estimates
upper_centile = floor(prob_coverage*size(sorted_data,1)); % upper bound
nCIs = size(sorted_data,1) - upper_centile;
HDI = zeros(2,grouping);

%% start
figure; hold on

%% select color scheme
color_scheme = cubehelixmap('semi_continuous',grouping+1);

%% Scatter plot of the data with automatic spread
% scatter plot parameters

for u=1:grouping
    tic
    tmp = sort(Data(~isnan(Data(:,u)),u));
    % create a matrix with spread = 2
    Y = NaN(length(tmp),2);
    Y(1,1) = tmp(1);
    c_index = [2 1];
    c_thresh = round(range(tmp)/50,1);
    for c = 2:length(tmp)
        if tmp(c)-tmp(c-1) < c_thresh
            Y(c,c_index(1)) = tmp(c);
            c_index = fliplr(c_index);
        else
            Y(c,1) = tmp(c);
            c_index = [2 1];
        end    
    end
    % find outliers
    class = rst_outlier(tmp,outlier_method);
    outliers = find(class);
    % plot
    X = repmat([0 within_gp_dispersion],[length(tmp),1]) + gp_index(u);
    X(isnan(Y)) = NaN;
    for p=1:size(Y,2)
        scatter(X(:,p),Y(:,p),point_size,color_scheme(u,:));
        scatter(X(outliers,p),Y(outliers,p),point_size,color_scheme(u,:),'filled');
    end
    toc
    %% Add the density estimate 
    % get the kernel
    [bc,K]=rst_RASH(tmp,100,dist_method);
    % remove 0s
    bc(K==0)=[]; K(K==0)= [];
    % create symmetric values
    K = (K - min(K)) ./ max(K); % normalize to [0 1] interval
    high=(K/2); low=(-high);
    
    % plot contours
    y1 = plot(high+gp_index(u),bc); set(y1,'Color',color_scheme(u,:)); hold on
    y2 = plot(low+gp_index(u),bc); set(y1,'Color',color_scheme(u,:));
    if isnumeric(y1)
        y1 = get(y1); y2 = get(y2); % old fashion matlab
    end
    % fill
    xpoints=[y2.XData',y1.XData']; filled=[y2.YData',y1.YData'];
    % check that we have continuous data, otherwise 0 pad
    if diff(xpoints(1,:)) > 0.001*(range(xpoints(:,1)))
        xtmp = NaN(size(xpoints,1)+1,size(xpoints,2));
        xtmp(1,:) = [gp_index(u) gp_index(u)];
        xtmp(2:end,:) = xpoints;
        xpoints = xtmp; clear xtmp
        ytmp = NaN(size(filled,1)+1,size(filled,2));
        ytmp(1,:) = filled(1,:);
        ytmp(2:end,:) = filled;
        filled = ytmp; clear ytmp        
    end
    
    if diff(xpoints(end,:)) > 0.001*(range(xpoints(:,1)))
        xpoints(end+1,:) = [gp_index(u) gp_index(u)];
        filled(end+1,:)  = filled(end,:);
    end   
    hold on; fillhandle=fill(xpoints,filled,color_scheme(u,:));
    set(fillhandle,'LineWidth',2,'EdgeColor',color_scheme(u,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color

    %% add IQR - using again Harell-Davis Q
    ql = rst_hd(tmp,0.25);
    [~,position] = min(abs(filled(:,1) - ql));
    plot(xpoints(position,:),filled(position,:),'Color',color_scheme(u,:));
    if  strcmpi(estimator,'median')
        qm = est(u);
    else
        qm = rst_hd(tmp,0.5);
    end
    [~,position] = min(abs(filled(:,1) - qm));
    plot(xpoints(position,:),filled(position,:),'Color',color_scheme(u,:));
    qu = rst_hd(tmp,0.75);
    [~,position] = min(abs(filled(:,1) - qu));
    plot(xpoints(position,:),filled(position,:),'Color',color_scheme(u,:));
    clear xpoints filled
    
    %% Hight Density Intervals
    tmp = sorted_data(:,u);
    ci = 1:nCIs; ciWidth = tmp(ci+upper_centile) - tmp(ci); % all centile distances
    [~,index]=find(ciWidth == min(ciWidth)); % densest centile
    if length(index) > 1; index = index(1); end % many similar values
    HDI(1,u) = tmp(index);
    HDI(2,u) = tmp(index+upper_centile);
    
    % plot this with a rectangle
    X = (gp_index(u)-0.3):0.1:(gp_index(u)+0.3);
    plot(X(3:5),repmat(est(u),[1,length(X(3:5))]),'LineWidth',6,'Color',color_scheme(u,:));
    rectangle('Position',[X(3),HDI(1,u),X(5)-X(3),HDI(2,u)-HDI(1,u)],'Curvature',[0.4 0.4],'LineWidth',3,'EdgeColor',color_scheme(u,:))
end


%% finish off
cst = max(abs(diff(Data(:)))) * 0.1;
axis([0.3 gp_index(grouping)+0.7 min(Data(:))-cst max(Data(:))+cst])    

if size(Data,1) == 1
    title(sprintf('Data distribution with %s and 95%% High Density Interval',estimator),'FontSize',16);
else
    title(sprintf('Data distributions with %ss and 95%% High Density Intervals',estimator),'FontSize',16);
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

if exist('plotly','file') == 2
    output = questdlg('Do you want to output this graph with Plotly?', 'Plotly option');
    if strcmp(output,'Yes')
        save_dir = uigetdir('select directory to save html files','save in');
        if ~isempty(save_dir)
            cd(save_dir)
            try
                fig2plotly(gcf,'strip',true,'offline', true);
            catch
                fig2plotly(gcf,'strip',true); % in case local lib not available
            end
        else
            return
        end
    else
        return
    end
end