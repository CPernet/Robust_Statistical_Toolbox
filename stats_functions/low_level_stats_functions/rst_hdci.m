function [xhd CI] = rst_hdci(x,nboot,plotCI)

% [xhd CI] = hdci(x,nboot,plotCI)
% HDCI computes the 95% confidence interval 
% for the deciles of a distribution
% using the Harrell-Davis estimate  
% See Wilcox p.130-134
%
% GAR, University of Glasgow, Sep 2009 

rand('state',sum(100*clock));

if nargin < 2;nboot=100;plotCI=0;end

n=length(x);                    

for d=1:9
    
    % The constant c was determined so that the
    % probability coverage of the confidence interval is
    % approximately 95% when sampling from normal and
    % non-normal distributions

    c = 1.96 + .5064.* (n.^-.25);

    if d<=2 || d>=8
        if n <=21
            c = -6.23./n+5.01;
        end
    end
    
    if d<=1 || d>=9
        if  n<=40
            c = 36.2./n+1.31;
        end
    end

if n<=10
    error('confidence intervals of the hd estimates of the deciles cannot be computed for less than 11 observations')
end
       
   xhd(d) = hd(x,d./10);
   xd_bse = bootse(x,nboot,'hd',d./10); % bse = bootse(x,nboot,est)
   CI(d,1) = xhd(d)-c.*xd_bse;
   CI(d,2) = xhd(d)+c.*xd_bse;
end

if plotCI==1
    figure;set(gcf,'Color','w');hold on
    plot(1:9,xhd,'k.',1:9,CI(:,1),'r+',1:9,CI(:,2),'r+')
%     for ref=1:9
%     refline(0,xhd(ref))
%     end

%     xlabel('x (first group)','FontSize',16)
%     ylabel('Delta','FontSize',16)
    set(gca,'FontSize',14)
    box on
end

return
