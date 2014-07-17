function thetaq = rst_hd(x,q)

% Compute the Harrell-Davis estimate of the qth quantile
% The vector x contains the data, and the desired quantile is q
% The default value for q is .5.
%
% Original R code by Rand Wilcox
% See Wilcox p.71, 139
% GAR, University of Glasgow, Dec 2007

if nargin<2; q=.5;end
n=length(x);
m1=(n+1).*q;
m2=(n+1).*(1-q);
vec=1:length(x);
w=betacdf(vec./n,m1,m2)-betacdf((vec-1)./n,m1,m2); 
        
y=sort(x);
thetaq=sum(w(:).*y(:));

return

% Data from Wilcox p.140
% heartbeat=[190 80 80 75 50 40 30 20 20 10 10 10 0 0 -10 -25 -30 -45 -60 -85];

% a=[77 87 88 114 151 210 219 246 253 262 296 299 306 376 428 515 666 1310 2611]; % Wilcox, Table 3.2 p.58
