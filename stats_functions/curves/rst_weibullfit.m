function fit = rst_weibullfit(AFC, stimlevel, pc, CritY, labelx, estdelta, flag)

% Fit a cumulative Weibull function to relationship between vectors:
% stimlevel and pc, both obtained in N-alternative discrimination task.
% Set AFC to N.
%
% CritY = threshold proportion correct. default is halfway between chance
% and 100-percent correct.
% labelx = string to label the x axis of automatically generated graph
% estdelta = see two paragraphs below.
%
% Set estdelta to [] to set delta to be a constant value of 0; e.g., assume
% that observer can achieve one-hundred percent accuracy at some
% arbitrarily high stimulus value.
%
% Otherwise delta is a free parameter (default), and estdelta is the initial
% estimate. Default value for estdelta (if not set) is accuracy for the
% highest stimulus level.
%
% ----- An example using mock data -----
%
% X=0:0.1:1; % mock stimulus levels
% p=[1 1 0.01]; % assumed "true" that parameters driving observer accuracy
% AFC=5;
% Gamma=1/AFC;
% wblpsy=inline('(1-g-p(3))*wblcdf(X,p(1),p(2))+g','p','g','X'); % need this to model performance
% Y=Scale(wblpsy(p,Gamma,X) + randn(1,length(X))*(1/100)); % actual performance is like model but with noise
% fit = weibullfit(AFC, X, Y, [.3 .6 .9], 'Nostril Expansion', 0.1 );
%
if ~AFC
    Gamma = 0;
else
    Gamma=1/AFC;
end

plotthreshpc=(1-Gamma)/2 + Gamma;
Alpha=1;
Beta=1;
[maxx,maxxloc]=max(stimlevel);
if nargin < 6    
    Delta=1-pc(maxxloc);
elseif nargin >= 6
    Delta=estdelta;
end



if ~isempty(Delta)
    Init=[Alpha Beta Delta];
    wblpsy=inline('(1-g-p(3))*wblcdf(X,p(1),p(2))+g','p','g','X');
    fiterr=inline('sum(( ((1-g-p(3))*wblcdf(X,p(1),p(2))+g) - Y ).^2)', 'p', 'g', 'X', 'Y');
    % Analytic solution for 2AFC
    %
    % threshold=inline('exp(log( (-p(1)^p(2))*log(1-((Y-0.5)/(0.5-p(3)))) )/p(2))','p','Y');
else
    Init=[Alpha Beta];
    wblpsy=inline('(1-g-0)*wblcdf(X,p(1),p(2))+g','p','g','X');
    fiterr=inline('sum(( ((1-g-0)*wblcdf(X,p(1),p(2))+g) - Y ).^2)', 'p', 'g', 'X', 'Y');    
end


% Obtain best-fitting parameters for psychometric function
stimlevel=stimlevel(:); pc=pc(:); CritY=CritY(:);
fit.params = fminsearch(@(p) fiterr(p,Gamma,stimlevel,pc), Init);
if ~isempty(Delta) && fit.params(3)<0
    % If delta is a free parameter, and
    % the estimated maximum for accuracy is above 1 then make delta a fixed
    % variable equal to zero.
    Init=[Alpha Beta];
    wblpsy=inline('(1-g-0)*wblcdf(X,p(1),p(2))+g','p','g','X');
    fiterr=inline('sum(( ((1-g-0)*wblcdf(X,p(1),p(2))+g) - Y ).^2)', 'p', 'g', 'X', 'Y');    
    fit.params = fminsearch(@(p) fiterr(p,Gamma,stimlevel,pc), Init);
    estdelta=[];
    Delta=[];
end

% Interpolate stimulus thresholds at criteria accuracy CritY
% InitLevel=mean(stimlevel);
[mv,mloc]=min(abs(pc-plotthreshpc));
InitLevel=stimlevel(mloc);
fit.criteriapc=CritY;
for i=1:length(fit.criteriapc)
    fit.thresholds(i) = fminsearch(@(thestimlevel) fiterr(fit.params,Gamma,thestimlevel,fit.criteriapc(i)), InitLevel);
end
 

% original plot code
if flag == 1
figure; hold on;
plot(stimlevel,pc,'ro','MarkerSize',12,'LineWidth',2); % plot data
X=linspace(min(stimlevel),maxx,100); % plot fitted function
plot(X,wblpsy(fit.params,Gamma,X),'k','MarkerSize',12,'LineWidth',2);
title(['Psychometric function for ',num2str(AFC),'-alternative discrimination'],'FontSize',14);
set(gca,'YLim',[Gamma 1]);
ylabel('Proportion correct','FontSize',14);
if nargin<5
    xlabel('Stimulus level X','FontSize',14);
else
    xlabel(labelx,'FontSize',14);
end

% plot delta for illustration
if ~isempty(Delta)
    plot([min(X) maxx],[1-fit.params(3) 1-fit.params(3)],'r--');
else
    plot([min(X) maxx],[1 1],'r--');
end
 
 
% Find and plot threshold at point of maximum inflection
threshold = fminsearch(@(thestimlevel) fiterr(fit.params,Gamma,thestimlevel,plotthreshpc), InitLevel);
plot([min(X) maxx],[plotthreshpc plotthreshpc],'b--');
if ~isempty(Delta)
    legend('data','fitted function',['delta is ',num2str(fit.params(3))],[num2str(round(plotthreshpc*100)),'-percent threshold at ',num2str(threshold),' units'],'Location','SouthEast');
else
    legend('data','fitted function','delta is 0',[num2str(round(plotthreshpc*100)),'-percent threshold at ',num2str(threshold),' units'],'Location','SouthEast');
end
plot([threshold threshold],[Gamma 1],'b--');
end

fit.predpc=wblpsy(fit.params,Gamma,stimlevel);


return