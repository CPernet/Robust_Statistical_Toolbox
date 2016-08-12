function [class,distance] = rst_outlier(data,method)

% Allows to detect outliers 
%
% FORMAT [class,distance] = rst_outlier(data,method)
%
% IMPUT data a vector of data points
%       method = 1 - is the standard MAD Median rule
%              = 2 - is the MAD Median rule with adustement for small samples
%              = 3 - is the S-outlier method
%
% OUTPUT class is a binary vector indicating outliers
%        distance is a metric reporting the relative distance of each data point
%
% Methods 1/2 - this is similar, in spirit, to reject data points located above 2
% standard deviations - however, using std is biased because by definition
% the mean and the std are affected by the outlier(s); instead one rely on
% the median and the deviation to the median (mad). We then apply a
% consistency factor so the mad estimates the std and find the outliers
% as data points above a given limit (the 75th quantile of the
% distribution)
%
% Method 3 allows to detect outliers in the data using a S estimator
% Like the MAD, the S estimator is robust and relies on the median of
% absolute distances. However, MAD relies on the distance to the median
% whereas S relies on the median of absolute distances. Because it doesn't
% rely on an estimator of central tendency (like the MAD) it works well for
% non symmetric distributons
%
% Ref. Rousseeuw, P.J. and Croux C. (1993). Alternatives to the the median
% absolute deviation. Journal of the American Statistical Association, 88
% (424) p 1273-1263
%
% Cyril Pernet v1 - May 2012
% -----------------------------

%% varargin
k = 2.2414; % = sqrt(chi2inv(0.975,1))
[n,p]=size(data);
if nargin < 2
    method = 2;
    disp('default method selected: MAD with adjustement')
elseif nargin > 2
    error('only two arguments in')
end

if ~isnumeric(method) ||  method > 3
    error(sprintf('method must be a numeric between 1 and 3: \n 1 MAD with finite sample correction \n 2 MAD \n 3 S-outliers'))
end

%% compute
if method == 1 || method == 2
    M = nanmedian(data);
    MAD=nanmedian(abs(data - repmat(nanmedian(data),n,1)));
end

switch method
    
    case 1 % Median Absolute Deviation with finite sample correction factor
        
        if n == 2
            bn=1.197;
        elseif n == 3
            bn=1.49;
        elseif n == 4
            bn=1.36;
        elseif n == 5
            bn=1.217;
        elseif n == 6
            bn=1.189;
        elseif n == 7
            bn=1.138;
        elseif n == 8
            bn=1.127;
        elseif n == 9
            bn=1.101;
        else
            bn=n/(n-0.8);
        end
        
        MADS=repmat((MAD.*1.4826.*bn),n,1);
        class = data > (repmat(M,n,1)+(k.*MADS));
        class = class+isnan(data);
        distance = MADS(1,:);
        
    case 2 % Normalized Median Absolute Deviation
        
        MADN = repmat((MAD./.6745),n,1); % same as MAD.*1.4826 :-)
        class = (abs(data-repmat(M,n,1)) ./ MADN) > k;
        class = class+isnan(data);
        distance = MADN(1,:);
        
    case 3 % S outliers
        
        distance = NaN(n,p);
        for p=1:size(data,2)
            tmp = data(:,p);
            points = find(~isnan(tmp));
            tmp(isnan(tmp)) = [];
            
            % compte all distances
            n = length(tmp);
            for i=1:n
                j = points(i);
                indices = [1:n]; indices(i) = [];
                distance(j,p) = median(abs(tmp(i) - tmp(indices)));
            end
            
            % get the S estimator
            % consistency factor c = 1.1926;
            Sn = 1.1926*median(distance(points,p));
            
            % get the outliers in a normal distribution
            class(:,p) = (distance(:,p) ./ Sn) > k; % no scaling needed as S estimates already std(data)
            class(:,p) = class(:,p)+isnan(data(:,p));
        end
end
