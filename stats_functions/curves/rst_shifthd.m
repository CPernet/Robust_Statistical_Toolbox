function [xd yd delta deltaCI] = rst_shifthd(x,y,nboot,plotshift)
% [xd yd delta deltaCI] = shifthd(x,y,nboot,plotshift)
% SHIFTHD computes the 95% simultaneous confidence intervals 
% for the difference between the deciles of two independent groups
% using the Harrell-Davis estimate  
% See Wilcox p.151-155
%
% GAR, University of Glasgow, Dec 2007 

rand('state',sum(100*clock));

if nargin < 3;nboot=200;plotshift=0;end

nx=length(x);
ny=length(y);
n=min(nx,ny);
c=(80.1./n.^2)+2.73; % The constant c was determined so that the simultaneous 
                     % probability coverage of all 9 differences is
                     % approximately 95% when sampling from normal
                     % distributions

for d=1:9
   xd(d) = hd(x,d./10);
   yd(d) = hd(y,d./10);
   xd_bse = bootse(x,nboot,'hd',d./10); % bse = bootse(x,nboot,est)
   yd_bse = bootse(y,nboot,'hd',d./10);
   delta(d) = yd(d)-xd(d);
   deltaCI(d,1) = yd(d)-xd(d)-c.*sqrt(xd_bse.^2+yd_bse.^2);
   deltaCI(d,2) = yd(d)-xd(d)+c.*sqrt(xd_bse.^2+yd_bse.^2);
end

if plotshift==1
    figure;set(gcf,'Color','w');hold on
    plot(xd,yd-xd,'k.',xd,deltaCI(:,1),'r+',xd,deltaCI(:,2),'r+')
    refline(0,0);
    xlabel('x (first group)','FontSize',16)
    ylabel('Delta','FontSize',16)
    set(gca,'FontSize',14)
    box on
end

% Data from Wilcox p.150
% control=[41 38.4 24.4 25.9 21.9 18.3 13.1 27.3 28.5 -16.9 26 17.4 21.8 15.4 27.4 19.2 22.4 17.7 26 29.4 21.4 26.6 22.7];
% ozone=[10.1 6.1 20.4 7.3 14.3 15.5 -9.9 6.8 28.2 17.9 -9 -12.9 14 6.6 12.1 15.7 39.9 -15.9 54.6 -14.7 44.1 -9];
