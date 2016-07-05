function [n,x,h] = rst_hist(data,varargin)

% rst_hist is a simple histogram plot with optimized number of bins by defailts and
% all the graphic options available (ie no need to set(get(gca), ... )
%
%  FORMAT [n,xout,h] = rst hist(data,bins, 'BarWidth',BarWidth, 'LineWidth',LineWidth, ...
%                     'FaceColor',FaceColor,'EdgeColor',EdgeColor,...
%                     'FaceAlpha',FaceAlpha, 'EdgeAlpha',EdgeAlpha)
% 
% INPUT data is a vector of data to plot
%       nbins is the number bins to use
%       'BarWidth' specifies the width of bars  
%       'FaceColor' specifies face color as RGB 
%       'EdgeColor' specifies edge color as as RGB 
%       'FaceAlpha' specifies transparency value of bars
%       'EdgeAlpha' specifies transparency value of bars' edges
%
% OUTPUT n is the frequency counts 
%        xout is the bin locations
%        h is the bin width
%
% Cyril Pernet August 2015
% ----------------------------------------------------------------------
% Copyright (C) RST Toolbox Team 2015

%% check arguments in

if ~isvector(data)
    error('data in must be a vector')
end

% set defaults
opt.BarWidth = 1;
opt.FaceColor = [0.5 0.5 1];
opt.EdgeColor = [0 0 0];
opt.FaceAlpha = 0.9;
opt.EdgeAlpha = 1;

if nargin>1
   for o=1:2:length(varargin)
       if strcmpi(cell2mat(varargin(o)),'BarWidth')
       opt.BarWidth = cell2mat(varargin(o+1));
       elseif strcmpi(cell2mat(varargin(o)),'FaceColor')
       opt.FaceColor = cell2mat(varargin(o+1));
       elseif strcmpi(cell2mat(varargin(o)),'EdgeColor')
       opt.EdgeColor = cell2mat(varargin(o+1));
       elseif strcmpi(cell2mat(varargin(o)),'FaceAlpha')
       opt.FaceAlpha = cell2mat(varargin(o+1));
       elseif strcmpi(cell2mat(varargin(o)),'EdgeAlpha')
       opt.EdgeAlpha = cell2mat(varargin(o+1));
       end
   end
end

%% basic stuff to make an histogram
data = data(~isnan(data));
k = round(1 + log2(length(data)));
[n,x]=hist(data,k);
h = x(2) - x(1);

%% plot
if nargout == 0
    try
        bar(x,n/(length(data)*h),opt.BarWidth, ...
        'FaceColor',opt.FaceColor,'EdgeColor',opt.EdgeColor, ...
        'FaceAlpha',opt.FaceAlpha,'EdgeAlpha',opt.EdgeAlpha);
    catch
        bar(x,n/(length(data)*h),opt.BarWidth, ...
        'FaceColor',opt.FaceColor,'EdgeColor',opt.EdgeColor);
    end
    xlabel('data values','FontSize',12); ylabel('Freq.','FontSize',12)
    grid on; hold on; set(gca,'Layer','Top')
end
