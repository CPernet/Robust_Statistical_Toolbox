function [xd yd delta deltaCI] = rst_shiftdhd(x,y,nboot,plotshift)
%[xd yd delta deltaCI] = shiftdhd(x,y,nboot,plotshift)
% shiftdhd computes the 95% simultaneous confidence intervals 
% for the difference between the deciles of two dependent groups
% using the Harrell-Davis estimate  
% See Wilcox p.184-187
%
% GAR, University of Glasgow, Dec 2007 

if nargin < 3;nboot=200;plotshift=0;end

n=length(x);
c=(37./n.^1.4)+2.75; % The constant c was determined so that the simultaneous 
                     % probability coverage of all 9 differences is
                     % approximately 95% when sampling from normal
                     % distributions
                     
% Get >>ONE<< set of B bootstrap samples
% The same set is used for all nine quantiles being compared
list = zeros(nboot,n);
for b=1:nboot
list(b,:) = randsample(1:n,n,true);
end

for d=1:9
   xd(d) = bat_hd(x,d./10);
   yd(d) = bat_hd(y,d./10);
   delta(d) = yd(d) - xd(d);
   for b=1:nboot
      bootdelta(b) = bat_hd(y(list(b,:)),d./10) - bat_hd(x(list(b,:)),d./10); 
   end
   delta_bse = std(bootdelta,0);
   deltaCI(d,1) = yd(d)-xd(d)-c.*delta_bse;
   deltaCI(d,2) = yd(d)-xd(d)+c.*delta_bse;
end

if plotshift==1
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
