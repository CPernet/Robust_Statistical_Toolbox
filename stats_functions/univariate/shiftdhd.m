function [xd, yd, delta, deltaCI] = shiftdhd(x,y,nboot,plotshift,list)
%[xd yd delta deltaCI] = shiftdhd(x,y,nboot,plotshift,list)
% shiftdhd computes the 95% simultaneous confidence intervals 
% for the difference between the deciles of two dependent groups
% using the Harrell-Davis estimate  
% See Wilcox p.184-187
%
% Original R code by Rand Wilcox
% GAR, University of Glasgow, Dec 2007 
% GAR, April 2013: replace randsample with randi
% GAR, June 2013: added list option

if nargin < 3;nboot=200;plotshift=0;list=[];end
if nargin < 4;plotshift=0;list=[];end
if nargin < 5;list=[];end

n=length(x);
c=(37./n.^1.4)+2.75; % The constant c was determined so that the simultaneous 
                     % probability coverage of all 9 differences is
                     % approximately 95% when sampling from normal
                     % distributions
                     
% Get >>ONE<< set of B bootstrap samples
% The same set is used for all nine quantiles being compared
if isempty(list)
    list = randi(n,nboot,n);
end

for d=1:9
   xd(d) = hd(x,d./10);
   yd(d) = hd(y,d./10);
   delta(d) = yd(d) - xd(d);
   for b=1:nboot
      bootdelta(b) = hd(y(list(b,:)),d./10) - hd(x(list(b,:)),d./10); 
   end
   delta_bse = std(bootdelta,0);
   deltaCI(d,1) = yd(d)-xd(d)-c.*delta_bse;
   deltaCI(d,2) = yd(d)-xd(d)+c.*delta_bse;
end

if plotshift==1
    figure;set(gcf,'Color','w');hold on
    plot(xd,yd-xd,'k.','MarkerSize',10)
    plot(xd,deltaCI(:,1),'r+',xd,deltaCI(:,2),'r+')
    refline(0,0);
    xlabel('x (first group)','FontSize',14)
    ylabel('Delta','FontSize',14)
    set(gca,'FontSize',12)
    box on
end

%% TEST
% Data from Wilcox p.187
% time1=[0 32 9 0 2 0 41 0 0 0 6 18 3 3 0 11 11 2 0 11];
% time3=[0 25 10 11 2 0 17 0 3 6 16 9 1 4 0 14 7 5 11 14];
