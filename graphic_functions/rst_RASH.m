function rst_RASH(data)

% compute the density estimate of data using RASH
% ie a Random Average Shifted Histogram algorithm
% Bourel et al. Computational Statistics and Data 
% Analysis 79 (2014) 149–164

n = length(data);
m = 100; % this is the parameter setting the number of hist in ASH and RASH

% standard probability histogram
% -----------------------------
% 1 determine how many bins
h = 2.15*sqrt(var(data))*n^(-1/5);
t0 = min(data)-1;
tm = max(data)+1;
bins = t0:h:tm;
% 2 compute frequency
nu=hist(data,bins);
nu(end) = [];
% 3 make it a prob
nu = nu./(n*h);
% 4 plot
bc = (t0+h/2):h:(tm-h/2);
figure; subplot(1,3,1);
bar(bc,nu,1,'FaceColor',[0.5 0.5 1]);
title('Prob histogram')

% general Average Shifted Histogram
% ----------------------------------
h = 2.15*sqrt(var(data))*n^(-1/5);
delta = h/m;
% 1 make a mesh with size delta
t0 = min(data)-1;
tf = max(data)+1;
nbin = ceil((tf-t0)/delta);
binedge = t0:delta:(t0+delta*nbin);
out = find(binedge>tf);
if out == 1
    binedge(out) = tf;
else
    binedge(out(1)) = tf;
    binedge(out(2:end)) = [];
end
% 2 bin count
nu = histc(data,binedge);
nu = [zeros(1,m-1) nu' zeros(1,m-1)];
% 3 Get the weight vector.
kern = inline('(15/16)*(1-x.^2).^2');
ind = (1-m):(m-1);
den = sum(kern(ind/m));% Get the denominator. 
wm = m*(kern(ind/m))/den;% Create the weight vector.
% 4 Get the bin heights over smaller bins.
ASH=zeros(1,nbin);
for k=1:nbin
   ind=k:(2*m+k-2);
   ASH(k)=sum(wm.*nu(ind));
end
ASH = ASH/(n*h);
bc = t0+((1:nbin)-0.5)*delta;
subplot(1,3,2); bar(bc,ASH,1,'FaceColor',[0.5 0.5 1]);
title(['ASH m = ' num2str(m)])

% RASH
% -----

h = 2.15*sqrt(var(data))*n^(-1/5);
delta = h/m;
% 1 make a mesh with size delta
t0 = min(data)-1;
tf = max(data)+1;
nbin = ceil((tf-t0)/delta);
binedge = t0:delta:(t0+delta*nbin);
if out == 1
    binedge(out) = tf;
else
    binedge(out(1)) = tf;
    binedge(out(2:end)) = [];
end
% 2 Get the weight vector.
kern = inline('(15/16)*(1-x.^2).^2');
ind = (1-m):(m-1);
den = sum(kern(ind/m));% Get the denominator. 
wm = m*(kern(ind/m))/den;% Create the weight vector.
% 3 compute bin with shited edges
RH=zeros(1,nbin);
RSH=zeros(m,nbin);
for e=1:m
    v = binedge + (delta*randn(1,1)); % e is taken from N(0,h);
    v(v<t0) = t0; % lower bound
    v(v>tf) = tf; % upper bound
    nu = histc(data,v);
    nu = [zeros(1,m-1) nu' zeros(1,m-1)];
    for k=1:nbin
        ind=k:(2*m+k-2);
        RH(k)=sum(wm.*nu(ind));
    end
    RSH(e,:) = RH/(n*h);
end

RASH = mean(RSH,1);
bc = t0+((1:nbin)-0.5)*delta;
subplot(1,3,3); bar(bc,RASH,1,'FaceColor',[0.5 0.5 1]);
title(['RASH m = ' num2str(m)])

