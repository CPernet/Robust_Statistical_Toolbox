function result = rst_Hotteling(Data,percent)


% test 
% N = 20;
% mu = [2 3 6];
% SIGMA = [1 0.5 0.2; 0.5 1 0.5; 0.2 0.5 1]; % data are correlated at .5 from 1 to 2 to 3 and .2 from 1 to 3
% for n=1:1000
% data = mvnrnd(mu,SIGMA,N);
% result = rst_Hotteling(data); % default 20% trimmed mean
% P(n) = result.p;
% end
% mean(P<0.05)

% one way Hotelling

if nargin < 2
    percent = 20/100;
end

% compute
[n,p]    = size(Data);
y        = rst_trimmean(Data,percent);   % robust means to compare
S        = rst_wincov(Data,percent);     % robut covariance to account for sphericity
C        = [eye(p-1) ones(p-1,1).*-1];   % contrast matrix
Tsquare  = n*(C*y')'*inv(C*S*C')*(C*y'); % usual Hotelling Tsquare

% get the stats
df   = p-1;                % number of levels -1
g    = 2*floor(percent*n); % number of items to trim
dfe  = (n-g)-p+1;          % usaully dfe = n-p+1 here use effective n

result.F     = (dfe/((n-g/2)*df)) * Tsquare; % usually ( dfe / ((n-1)*df) ) * Tsquare; 
result.p     = 1 - fcdf(result.F, df, dfe);
