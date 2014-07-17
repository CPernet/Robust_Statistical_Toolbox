function distribution = rst_vincentize(data,quantile,flag)

% compute and return a vincentized distribution

sorted_data = sort(data);
nbin = floor(length(data) / quantile); 
index1 = 1; 
index2 = nbin;

for i=1:(quantile-1)
    distribution{i} = sorted_data(index1:index2);  % data per quantile
    index1 = index2+1; index2 = index2+nbin;
end
distribution{quantile} = sorted_data(index1:end); % outside the loop, might be a few value larger than [index1:index2]

if flag ~= 0 
    v=[]; t=[];
    for i=1:quantile
        v=[v;mean(distribution{i})];
        t=[t;distribution{i}];
    end
    figure('Name','Vincentized data');
    subplot(1,2,1); bar(v); title('cumulative plot','Fontsize',14); 
    axis tight; grid on; ylabel('data value'); xlabel('bins')
    subplot(1,2,2); hist(t,quantile); title('freq per quantile','Fontsize',14); 
    axis tight; grid on; xlabel('data value'); ylabel('freq')
end
 