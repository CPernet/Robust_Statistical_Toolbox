function [xd,yd,delta,deltaCI] = rst_shifthd(varargin)

% Shift function analysis
%
% FORMAT: [xd,yd,delta,deltaCI] = shiftdhd(x,y,'group',value,'nboot',value,'plotshift', 'yes/no')
%
% INPUTS x and y are vectors - the distributions to compare
%        group is either dependent or independent
%        nboot is the number of bootraps to perform (200 default)
%        plotshit to plot the results of not ('yes' as default)
%
% OUTPUTS xd and yd are the Harrell-Davis estimates of deciles
%         delta the difference between xd and yd
%         deltaCI the 95% simultaneous confidence intervals for the difference
%
% GAR, University of Glasgow, Dec 2007
% Cyril Pernet - 2020 RST toolbox cleanup and optimization

%% inputs
if nargin < 4
    error('not enough input arguments')
else
    x = varargin{1};
    y = varargin{2};
end

if nargin < 5
    nboot     = 200;
    plotshift = 'yes';
else
    for v=3:nargin
        if ischar(varargin{v})
            if contains(varargin{v},'boot','IgnoreCase',true)
                nboot = varargin{v+1};
            elseif contains(varargin{v},'plot','IgnoreCase',true)
                plotshift = varargin{v+1};
            elseif contains(varargin{v},'group','IgnoreCase',true)
                group = varargin{v+1};
            end
        end
    end
end

%% compute
rng('shuffle')
% The constant c was determined so that the simultaneous
% probability coverage of all 9 differences is
% approximately 95% when sampling from normal
% distributions
    
if contains(group,'dependent','IgnoreCase',true)
    n = length(x);
    c = (37./n.^1.4)+2.75; 
    
    % The same set is used for all nine quantiles being compared
    list = zeros(nboot,n);
    for b=1:nboot
        list(b,:) = randsample(1:n,n,true);
    end
    
    for d=9:-1:1
        xd(d)        = rst_hd(x,d./10);
        yd(d)        = rst_hd(y,d./10);
        delta(d)     = yd(d) - xd(d);
        bootdelta    = rst_hd(y(list),d/10) - rst_hd(x(list),d/10);
        delta_bse    = std(bootdelta,0);
        deltaCI(d,1) = yd(d)-xd(d)-c.*delta_bse;
        deltaCI(d,2) = yd(d)-xd(d)+c.*delta_bse;
    end
    
elseif contains(group,'independent','IgnoreCase',true)
    nx = length(x);
    ny = length(y);
    n  = min(nx,ny);
    c  = (80.1./n.^2)+2.73; 
    
    for d=9:-1:1
        xd(d)        = rst_hd(x,d./10);
        yd(d)        = rst_hd(y,d./10);
        xd_bse       = rst_bootse(x,'harrell-davis',nboot,d/10); % bse = bootse(x,nboot,est)
        yd_bse       = rst_bootse(y,'harrell-davis',nboot,d/10);
        delta(d)     = yd(d)-xd(d);
        deltaCI(d,1) = yd(d)-xd(d)-c.*sqrt(xd_bse.^2+yd_bse.^2);
        deltaCI(d,2) = yd(d)-xd(d)+c.*sqrt(xd_bse.^2+yd_bse.^2);
    end
    
else
    error('''group'' argument must be ''dependent'' or ''independent''')
end

%% figure
if strcmpi(plotshift,'Yes')
    figure;set(gcf,'Color','w');
    
    subplot(2,1,1); hold on
    [bc,K]=rst_RASH(x,100,'ASH');
    bar(bc,K,1,'FaceColor',[1 0 0],'EdgeColor',[0 0 0],'FaceAlpha',0.5,'EdgeAlpha',0);
    [bc,K]=rst_RASH(y,100,'ASH');
    bar(bc,K,1,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.5,'EdgeAlpha',0);
    axis tight; grid on; box on; title('Smoothed Histograms')
    ylabel('freqency'); xlabel('observations')
    
    subplot(2,1,2); hold on
    plot(xd,yd-xd,'k','LineWidth',3)
    plot(xd,deltaCI(:,1),'k+',xd,deltaCI(:,2),'k+','LineWidth',2)
    hold on; fillhandle=fill([xd fliplr(xd)],[deltaCI(:,1)' fliplr(deltaCI(:,2)')],[0 0 0]);
    set(fillhandle,'LineWidth',2,'EdgeColor',[0 0 0],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    h = refline(0,0); set(h,'LineStyle','--','Color',[0 0 0],'LineWidth',3);
    xlabel('gp1 deciles (arbitrary unit)','FontSize',12)
    ylabel('gp1 - gp2 ','FontSize',12)
    grid on; box on
end


