%% *Testing of robust repeated measure ANOVA via Hotelling T2^
% We start by simulation under H0, taking normal, uncorrelated data and
% test a standard linear model, the standard Hotelling test and the robust
% version. 

nmc        = 10000;
p_robust   = NaN(nmc,1);
p_standard = NaN(nmc,1);
pb_robust   = NaN(nmc,1);
pb_standard = NaN(nmc,1);

% for ech monte carlo collect the p values
for MC = 1:nmc 
    fprintf('Running MC %g \n',MC);
    Data = randn(30,4); 
    RT2 = rst_rep_anova_T2(Data,ones(size(Data,1),1),4,1000,'',[],0);
    p_standard(MC) = RT2.p; 
    pb_standard(MC) = RT2.pboot; 
    RT2R = rst_rep_anova_T2(Data,ones(size(Data,1),1),4,1000,'',[],1);
    p_robust(MC) = RT2.p;
    pb_robust(MC) = RT2.pboot;
end

% compute the type 1 error, ie how many false positives
fprintf('the standard hotelling test for alpha=5%% gives %g%%\n',sum(p_standard<0.05)/MC*100)
fprintf('the robust hotelling test for alpha=5%% gives %g%%\n',sum(p_robust<0.05)/MC*100)

fprintf('the standard hotelling test for alpha=5%% gives pboot %g%%\n',sum(pb_standard<0.05)/MC*100)
fprintf('the robust hotelling test for alpha=5%% gives pboot %g%%\n',sum(pb_robust<0.05)/MC*100)

%% simply redo the dame adding a correlation between measurements (non sphericity) and outliers
nmc        = 10000;
p_robust   = NaN(nmc,1);
p_standard = NaN(nmc,1);
pb_robust   = NaN(nmc,1);
pb_standard = NaN(nmc,1);

% for ech monte carlo collect the p values
parfor MC = 1:nmc 
    fprintf('Running MC %g \n',MC); 
    Data = randn(30,4); 
    Data(1,:) = sort(Data(1,:)-2);  
    Data(2,:) = sort(Data(2,:))-2;
    Data(30,[3 4]) = Data(30,[3 4]) + [2 3];
    Data(29,[3 4]) = Data(29,[3 4]) + [2 3];
    RT2 = rst_rep_anova_T2(Data,ones(size(Data,1),1),4,1000,'',[],0);
    p_standard(MC) = RT2.p; 
    pb_standard(MC) = RT2.pboot; 
    RT2R = rst_rep_anova_T2(Data,ones(size(Data,1),1),4,1000,'',[],1);
    p_robust(MC) = RT2.p;
    pb_robust(MC) = RT2.pboot;
end

% compute the type 1 error, ie how many false positives
fprintf('the standard hotelling test for alpha=5%% gives %g%%\n',sum(p_standard<0.05)/MC*100)
fprintf('the robust hotelling test for alpha=5%% gives %g%%\n',sum(p_robust<0.05)/MC*100)

fprintf('the standard hotelling test for alpha=5%% gives pboot %g%%\n',sum(pb_standard<0.05)/MC*100)
fprintf('the robust hotelling test for alpha=5%% gives pboot %g%%\n',sum(pb_robust<0.05)/MC*100)

%% Redo the same with a 2 x 2
nmc = 1000
p_boot     = NaN(nmc,3);
p_robust   = NaN(nmc,3);
p_standard = NaN(nmc,3);

% for ech monte carlo collect the p values
for MC = 1:nmc 
    fprintf('Running MC %g \n',MC); 
    Data = randn(30,4); 
    Data(30,[2 4]) = Data(30,[2 4]) + [2 3];
    Data(29,[2 4]) = Data(29,[2 4]) + [2 3];
    RT2 = rst_rep_anova_T2(Data,ones(size(Data,1),1),[2 2],1000,'',[],0);
    p_standard(MC,:) = RT2.p; p_boot(MC,:) = RT2.pboot;
    RT2R = rst_rep_anova_T2(Data,ones(size(Data,1),1),[2 2],1000,'',[],1);
    p_robust(MC,:) = RT2.pboot;
end

% compute the type 1 error, ie how many false positives
p_standard=sum(p_standard<0.05)/nmc*100;
p_boot = sum(p_boot<0.05)/nmc*100;
p_robust = sum(p_robust<0.05)/nmc*100;
fprintf('the standard hotelling test for alpha=5%% gives for each  term %g%% %g%% %g%%\n',p_standard(1),p_standard(2),p_standard(3))
fprintf('the standard hotelling test with bootstrap, for alpha=5%% gives %g%% %g%% %g%%\n',p_boot(1),p_boot(2),p_boot(3))
fprintf('the robust hotelling test with bootstrap, for alpha=5%% gives %g%% %g%% %g%%\n',p_robust(1),p_robust(2),p_robust(3))






