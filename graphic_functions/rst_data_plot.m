function [est,HDI]=rst_data_plot(Data,varargin)

% plots the data split by groups showing each data points with the
% distribution and a summary statistics estimator with 95&% Bayes boot HDI
%
% FORMAT: [est,HDI]= rst_data_plot(Data,options)
%         [est,HDI]= rst_data+plot(Data,'between',0.25,'within',0.025,'pointsize',50,'estimator','median') 
%
% INPUT: Data is a matrix, data are plotted colmun-wise
%        options are 
%                'between' for the distance between distributions
%                'within' for the distance/width between points in a group
%                'point_size' the size of data points
%                'estimator' can be 'median', 'mean', 'trimmed mean' 
% OUTPUT: est is the summary statistics
%         HDI the 95% High Density Interval (Bayes bootstrap)
%
% example: Data = [randn(100,1) randg(1,100,1) randn(100,1) 1-randg(1,100,1)];
%          Data(50:55,1) = NaN; Data(90:100,2) = NaN;
%          [est,HDI]=rst_data_plot(Data,'between',0.25,'within',0.025,'pointsize',50,'estimator','median')
%
% see also cubehelixmap, rst_outlier, rst_RASH, rst_trimmean, rst_hd
%
% Cyril Pernet - the university of Edinburgh
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
estimator = 'HD';

% check inputs
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

%% start
figure; hold on

%% how many groups
grouping = size(Data,2);
gp_index = 1:(1+between_gp_dispersion):(grouping*(1+between_gp_dispersion));

%% select color scheme
color_scheme = cubehelixmap('semi_continuous',grouping+1);

%% Scatter plot of the data with automatic spread
% scatter plot parameters

for u=1:grouping
    tmp = sort(Data(~isnan(Data(:,u)),u));
    % find a spread needed
    s = sum(diff(tmp) < 0.1);
    spread = [0:(1/s):1];
    % creater a matrix with that spread
    change = find(diff(tmp) < 0.1);
    Y = NaN(length(tmp),2);
    Y(1:change(1),1) = tmp(1:change(1));
    c_index = 2;
    for c=2:length(change)
        Y((change(c-1)+1):change(c),c_index) = tmp((change(c-1)+1):change(c));
        if mod(c,2) == 0
            c_index = 1;
        else
            c_index = 2;
        end
    end
    Y((change(c)+1):length(tmp),c_index) = tmp((change(c)+1):length(tmp));
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
end   
cst = max(abs(diff(Data(:)))) * 0.1;
axis([0.3 gp_index(grouping)+0.7 min(Data(:))-cst max(Data(:))+cst])    
    
%% Add the density estimate 
for u=1:grouping
    tmp = sort(Data(~isnan(Data(:,u)),u));
    % get the kernel
    [bc,K]=rst_RASH(tmp,100,dist_method);
    % remove 0s
    bc(find(K==0)) = [];
    K(find(K==0))= [];
    % map bc onto observed values
    steps = range(tmp)/(length(bc)-1);
    newbc = min(tmp):steps:max(tmp);
    % create symmetric values 
    K = K/max(K); high=(K/2); low=(-high);
    % plot
    y1 = plot(high+gp_index(u),newbc); set(y1,'Color',color_scheme(u,:)); hold on
    y2 = plot(low+gp_index(u),newbc); set(y1,'Color',color_scheme(u,:)); 
    xpoints=[y2.XData',y1.XData']; filled=[y2.YData',y1.YData'];
    hold on; fillhandle=fill(xpoints,filled,color_scheme(u,:));
    set(fillhandle,'EdgeColor',color_scheme(u,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color   
end

%% Add the summary stat with 95% HDI (Bayes bootstrap)
if strcmpi(estimator,'Mean')
    est = nanmean(Data,1);
elseif strcmpi(estimator,'Trimmed mean')
    est = rst_tm(Y,trimming); 
elseif strcmpi(estimator,'Median')
    est = rst_hd(Data,decile); % Median estimated using the 5th decile of Harell Davis
end


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
        bb(boot,:) = rst_tm(resample,20); 
    elseif strcmpi(estimator,'Median')
        bb(boot,:) = rst_hd(resample,.5); 
    end
end

sorted_data = sort(bb,1); % sort bootstrap estimates
upper_centile = floor(prob_coverage*size(sorted_data,1)); % upper bound
nCIs = size(sorted_data,1) - upper_centile;
HDI = zeros(2,grouping);

for u = 1:grouping
    tmp = sorted_data(:,u);
    ci = 1:nCIs; ciWidth = tmp(ci+upper_centile) - tmp(ci); % all centile distances
    [~,index]=find(ciWidth == min(ciWidth)); % densest centile
    if length(index) > 1; index = index(1); end % many similar values
    HDI(1,u) = tmp(index);
    HDI(2,u) = tmp(index+upper_centile);
    
    % plot this with a rectangle
    X = [(gp_index(u)-0.3):0.1:(gp_index(u)+0.3)];
    plot(X(3:5),repmat(est(u),[1,length(X(3:5))]),'LineWidth',6,'Color',color_scheme(u,:));
    rectangle('Position',[X(3),HDI(1,u),X(5)-X(3),HDI(2,u)-HDI(1,u)],'Curvature',[0.4 0.4],'LineWidth',3,'EdgeColor',color_scheme(u,:))
end


%% finish off
if size(Data,1) == 1
title(sprintf('Data distribution with %s and 95%% High Dentity Interval',estimator));
else
   title(sprintf('Data distributions with %ss and 95%% High Dentity Intervals',estimator));
end
grid on
box on

