function[x,nu,xp,yp] = rst_density_hist(data,units)

% simply function returning the density histogram
% function adapted from the computational statistics toolbox
% http://www.pi-sigma.info
% 
% INPUT = data  - vector of data to plot
%       = 'units' - label for the x axis
% OUPUT x = values of each bin
%       nu = nb of observation per bin
%       xp = values for each bin of the pdf
%       yp = nb of expected observations per bin
%
% Cyril Pernet 18-Feb-2011
% ----------------------------------------------------------------------
% Copyright (C) RST Toolbox Team 2015

if nargin == 1
    units = 'data values';
end

% get mean and var
mu = nanmean(data);
v = nanvar(data);

% get the normal pdf for this distribution
xp = linspace(min(data),max(data));
if v <= 0
   error('Variance must be greater than zero')
   return
end
arg = ((xp-mu).^2)/(2*v);
cons = sqrt(2*pi)*sqrt(v);
yp = (1/cons)*exp(-arg);

% get histogram info 
k = round(1 + log2(length(data)));
[nu,x]=hist(data,k);
h = x(2) - x(1);

% plot
if nargout == 0
    % figure; set(gcf,'Color','w');
    bar(x,nu/(length(data)*h),1,'FaceColor',[0.5 0.5 1]);
    grid on; hold on; % axis('tight');
    plot(xp,yp,'r','LineWidth',3);
    xlabel(units,'FontSize',10); ylabel('Freq.','FontSize',10)
    title({['Density Histogram']; ['and Density Estimate']},'FontSize',12); % ,'FontWeight','demi')
end


