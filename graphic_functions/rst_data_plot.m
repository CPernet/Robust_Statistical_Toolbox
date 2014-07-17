function rst_data_plot(data,gp)

% plots the data split by levels and groups
% data is a matrix such as
%          size(data,1) = length(gp)
%          size(data,2) = levels
% for instance if size(data) = [20,2] we can define a vector gp and the
% plot will be data split per gp and in each gp 2 subplots
% there is one boxplot with median and interquartile range and there is a
% bar graph with the within-subject variance bar in each group - the same
% is plotted as curves

% check input
% ------------
if nargin >2 
    error('wrong number of arguments');
elseif nargin == 1
    N = size(data,1);
    levels = size(data,2);
    gp = ones(N,1); v = 1; n = N;
else
    N = size(data,1);
    levels = size(data,2);
    v = unique(gp);
    for i=1:length(v)
        n(i) = sum(gp == v(i));
    end
end

% histograms
% ----------
figure('Name','histrograms');
index = 1;
for i=1:length(v)
    tmp = data(gp == v(i),:);
    for j=1:levels
        subplot(length(v),levels,index);
        hist(tmp(:,j)); grid on
        title(['data gp ' num2str(i) 'condition' num2str(j)],'FontSize',14);
     index = index + 1;
   end
end


% boxplot
% ---------
rows = ceil(length(v)/3); % max 3 subplots per rows
columns = ceil(length(v) / rows); 
figure('Name','boxplots');
for i=1:length(v)
    tmp = data(gp == v(i),:);
    subplot(rows,columns,i);
    boxplot(tmp,'outliersize', 10); grid on
    title(['data gp ' num2str(i)],'FontSize',14);
end

% bar graphs
% ----------
figure('Name','Bars and standard errors');
for i=1:length(v)
    tmp = data(gp == v(i),:);
    subplot(rows,columns,i);
    stderror = std(tmp)/sqrt(n(i));
    bar(mean(tmp)); hold on
    errorbar(nanmean(tmp),stderror,'r','LineWidth',2); 
    title(['data gp ' num2str(i)],'FontSize',14);grid on
end

% curves
% ---------
figure('Name','Bars and standard errors');
for i=1:length(v)
    tmp = data(gp == v(i),:);
    subplot(rows,columns,i);
    stderror = std(tmp)/sqrt(n(i));
    plot(mean(tmp)); hold on
    errorbar(nanmean(tmp),stderror,'r','LineWidth',2); 
    title(['data gp ' num2str(i)],'FontSize',14);grid on
end

