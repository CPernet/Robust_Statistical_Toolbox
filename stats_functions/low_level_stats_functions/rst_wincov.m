function wincov = rst_#wincov(Data,Percent)

% Computes the winsorized covariance between the columns of the matrix Data
%
% FORMAT
%
% INPUT
%
% OUTPUT
%
% Cyril Pernet June 2012

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
low_outliers = NaN(g+1,size(Data,2));
high_outliers = NaN(g+1,size(Data,2));

for i=1:size(Data,2)
    low_outliers(:,i) = Sorted_Data(find(Sorted_Data(:,i)<=repmat(lower_bounds(i),n,1)),i);
    high_outliers(:,i) = Sorted_Data(find(Sorted_Data(:,i)>=repmat(higher_bounds(i),n,1)),i);
end

% substitute in original data to keep the relationship and compute covariance
WinData = Data;
for i=1:size(Data,2)
    for j=1:(g+1)
        WinData(Data(:,i) == low_outliers(j,i),i) = lower_bounds(i);
        WinData(Data(:,i) == high_outliers(j,i),i) = higher_bounds(i);
    end
end

% % we can check it works
% Sorted_WinData=NaN(size(Data));
% for i=1:size(Data,2)
%     Sorted_WinData(:,i) = sort(squeeze(WinData(:,i)));
% end
% [Sorted_Data Sorted_WinData] 
 
wincov = cov(WinData);

