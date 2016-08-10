function rst_data_plot(data)

% plots the data split by groups

Data = randn(100,3);

figure; hold on

%% how many groups
grouping = size(Data,2);

%% select color scheme
% intensity based: Green, D. A., 2011, `A colour scheme for the display of
% astronomical intensity images', Bulletin of the Astronomical Society of 
% India, 39, 289.

color_scheme = cubehelix(grouping+1,[0.5,-1.5,1,1], [0,1], [0,1]);

%% Scatter plot of the data with automatic spread
% scatter plot parameters
within_gp_dispersion = 0.025;
point_size = 50;

for u=1:grouping
    tmp = sort(Data(~isnan(Data(:,u)),u));
    % find a spread needed
    s = sum(diff(tmp) < 0.1);
    spread = [0:(1/s):1];
    % creater a matrix with that spread
    change = find(diff(tmp) < 0.1);
    Y = NaN(length(tmp),2);
    Y(1:change(1),1) = tmp(1:change(1));
    c_index = 2;
    for c=2:length(change)
        Y((change(c-1)+1):change(c),c_index) = tmp((change(c-1)+1):change(c));
        if mod(c,2) == 0
            c_index = 1;
        else
            c_index = 2;
        end
    end
    Y((change(c)+1):length(tmp),c_index) = tmp((change(c)+1):length(tmp));
    % plot
    X = repmat([0 within_gp_dispersion],[length(tmp),1]) + u;
    X(isnan(Y)) = NaN;
    for p=1:size(Y,2)
        scatter(X(:,p),Y(:,p),point_size,color_scheme(u,:));
    end
end   
cst = max(abs(diff(Data(:)))) * 0.1;
axis([0.5 grouping+0.5 min(Data(:))-cst max(Data(:))+cst])    
    
%% Add the density estimate 
for u=1:grouping
    tmp = sort(Data(~isnan(Data(:,u)),u));
    [N,X]=rst_RASH(tmp);
    % plot
    
    
end

%% Add the summary stat with 95% HDI (Bayes bootstrap)


%% finish off
grid on
box on

