function out = rst_dynamicdprime(data,quantile)

% compute the dynamic d prime based on the vincentization of RT
% data is a matrix of RT (column 1) hits (column2) and stimulus codes (column3)
% quantile is the number of bin to do
% out is the result

% take stimulus info
% ---------------------
Stim = data(:,3); 

% sort RT per stimulus
% -------------------
% apply a filter to find relevant stimuli - here odd/even numbers
stim1 = find(mod(Stim,2) == 1);
stim2 = find(mod(Stim,2) == 0);
% sort RT and get there index
[stim1RT,stim1index] = sort(data(stim1,1));
[stim2RT,stim2index] = sort(data(stim2,1));
% apply the RT sorting to Perf
stim1Perf = data(stim1,2);
stim1Perf = stim1Perf(stim1index);
stim2Perf = data(stim2,2);
stim2Perf = stim2Perf(stim2index);

% get the nb of items per bin per stimulus
% -----------------------------------------
nbin1 = floor(length(stim1RT) / quantile); 
index1 = 1; index2 = nbin1;
nbin2 = floor(length(stim2RT) / quantile); 
index3 = 1; index4 = nbin2;

for i=1:(quantile-1)
    v1{i} = stim1Perf(index1:index2);  
    index1 = index2+1; index2 = index2+nbin1;
    v2{i} = stim2Perf(index3:index4);  
    index3 = index4+1; index4 = index4+nbin2;
end
v1{quantile} = stim1Perf(index1:end); % outside the loop, might be a few value larger than [index1:index2]
v2{quantile} = stim2Perf(index3:end);

% compute d' for each bin
% -------------------------
for i=1:quantile
    [zHits,zFA,out(i),c,cprime,beta] = rst_sdt1(length(v1{i}),length(v2{i}),sum(v1{1}),length(v2{i})-sum(v2{i}),0);
end



