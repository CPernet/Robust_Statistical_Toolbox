function wincov = rst_wincov(Data,Percent)

% Computes the winsorized covariance between the columns of the matrix Data
%
% FORMAT wincov = rst_wincov(Data,Percent)
%
% INPUT Data is a vector
%       Percent is the anount of trimming (default is 20%)
%
% OUTPUT wincov is the winsorized covariance matrix (ie the covariance of
%        the wonsorzed data)
%
% Cyril Pernet 
% -------------------------------------------------------------------------
% Copyright (C) RST Toolbox Team 2015

if nargin == 1
    Percent = 20/100;
end

if Percent > 1
    Percent = Percent / 100;
end

n = length(Data);
g=floor(Percent*n);

% sort data per column
Sorted_Data=NaN(size(Data));
for i=1:size(Data,2)
    Sorted_Data(:,i) = sort(squeeze(Data(:,i)));
end

lower_bounds=Sorted_Data(g+1,:);
higher_bounds=Sorted_Data(n-g,:);

% find which values have to be changed
for i=1:size(Data,2)
    low_outliers{i}  = Sorted_Data(find(Sorted_Data(:,i)<=repmat(lower_bounds(i),n,1)),i);
    high_outliers{i} = Sorted_Data(find(Sorted_Data(:,i)>=repmat(higher_bounds(i),n,1)),i);
end

% substitute in original data to keep the relationship and compute covariance
WinData = Data;
for i=1:size(Data,2)
    for j=1:length(low_outliers{i})
        WinData(Data(:,i) == low_outliers{i}(j),i) = lower_bounds(i);
    end
    
    for j=1:length(high_outliers{i})
        WinData(Data(:,i) == high_outliers{i}(j),i) = higher_bounds(i);
    end
end

wincov = cov(WinData);

