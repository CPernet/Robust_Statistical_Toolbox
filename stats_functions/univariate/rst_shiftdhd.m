function [xd,yd,delta,deltaCI] = rst_shiftdhd(varargin)

% Shift function analysis for dependent groups
%
% FORMAT: [xd,yd,delta,deltaCI] = shiftdhd(x,y,'nboot',value,'plotshift', 'yes/no')
%
% INPUTS x and y are vectors - the distributions to compare
%        nboot is the number of bootraps to perform (200 default)
%        plotshit to plot the results of not ('yes' as default)
%
% OUTPUTS xd and yd are the Harrell-Davis estimates of deciles
%         delta the difference between xd and yd
%         deltaCI the 95% simultaneous confidence intervals for the difference
%
% GAR, University of Glasgow, Dec 2007 
% Cyril Pernet - 2020 RST toolbox cleanup and optimnization

%% inputs
if nargin < 2
    error('not enough input arguments')
else
    x = varargin{1};
    y = varargin{2};
end

if nargin < 3
    nboot     = 200;
    plotshift = 'yes';
else
    for v=3:nargin
        if ischar(varargin{v})
            if contains(varargin{v},'boot','IgnoreCase',true)
                nboot = varargin{v+1};
            elseif contains(varargin{v},'plot','IgnoreCase',true)
                plotshift = varargin{v+1};
            end
        end
    end
end

%% compute

n=length(x);
c=(37./n.^1.4)+2.75; % The constant c was determined so that the simultaneous 
                     % probability coverage of all 9 differences is
                     % approximately 95% when sampling from normal
                     % distributions
                     
% The same set is used for all nine quantiles being compared
list = zeros(nboot,n);
for b=1:nboot
    list(b,:) = randsample(1:n,n,true);
end

for d=9:-1:1
    xd(d)    = rst_hd(x,d./10);
    yd(d)    = rst_hd(y,d./10);
    delta(d) = yd(d) - xd(d);
    parfor b=1:nboot
        bootdelta(b) = rst_hd(y(list(b,:)),d./10) - rst_hd(x(list(b,:)),d./10);
    end
    delta_bse = std(bootdelta,0);
    deltaCI(d,1) = yd(d)-xd(d)-c.*delta_bse;
    deltaCI(d,2) = yd(d)-xd(d)+c.*delta_bse;
end

if strcmpi(plotshift,'Yes')
    figure;set(gcf,'Color','w');hold on
    plot(xd,yd-xd,'k.',xd,deltaCI(:,1),'r+',xd,deltaCI(:,2),'r+')
    refline(0,0);
    xlabel('x (first group)','FontSize',14)
    ylabel('Delta','FontSize',14)
    set(gca,'FontSize',12)
    box on
end

% Data from Wilcox p.187
% time1=[0 32 9 0 2 0 41 0 0 0 6 18 3 3 0 11 11 2 0 11];
% time3=[0 25 10 11 2 0 17 0 3 6 16 9 1 4 0 14 7 5 11 14];
