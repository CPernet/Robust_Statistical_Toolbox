function [xd yd delta deltaCI] = rst_shifthd(x,y,nboot,plotshift)
% [xd yd delta deltaCI] = shifthd(x,y,nboot,plotshift)
% SHIFTHD computes the 95% simultaneous confidence intervals 
% for the difference between the deciles of two independent groups
% using the Harrell-Davis estimate  
% See Wilcox p.151-155
%
% GAR, University of Glasgow, Dec 2007 

rng('shuffle')
if nargin < 3;nboot=200;plotshift=1;end

nx=length(x);
ny=length(y);
n=min(nx,ny);
c=(80.1./n.^2)+2.73; % The constant c was determined so that the simultaneous 
                     % probability coverage of all 9 differences is
                     % approximately 95% when sampling from normal
                     % distributions

for d=1:9
   xd(d) = rst_hd(x,d./10);
   yd(d) = rst_hd(y,d./10);
   xd_bse = rst_bootse(x,'harrell-davis',nboot,d./10); % bse = bootse(x,nboot,est)
   yd_bse = rst_bootse(y,'harrell-davis',nboot,d./10);
   delta(d) = yd(d)-xd(d);
   deltaCI(d,1) = yd(d)-xd(d)-c.*sqrt(xd_bse.^2+yd_bse.^2);
   deltaCI(d,2) = yd(d)-xd(d)+c.*sqrt(xd_bse.^2+yd_bse.^2);
end

if plotshift==1
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

% Data from Wilcox p.150
% control=[41 38.4 24.4 25.9 21.9 18.3 13.1 27.3 28.5 -16.9 26 17.4 21.8 15.4 27.4 19.2 22.4 17.7 26 29.4 21.4 26.6 22.7];
% ozone=[10.1 6.1 20.4 7.3 14.3 15.5 -9.9 6.8 28.2 17.9 -9 -12.9 14 6.6 12.1 15.7 39.9 -15.9 54.6 -14.7 44.1 -9];
