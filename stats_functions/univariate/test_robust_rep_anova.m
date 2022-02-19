function test_robust_rep_anova


%% One sample repeated measure
% ----------------------------

N             = 20;                                        % 20 subjects
gp            = ones(N,1);                                 % 1 group
factors       = 3;                                         % N repeated meaures

% test under H0, to check the type 1 error rate
nboot = 10000;
P1    = NaN(nboot,1);
P2    = NaN(nboot,1);
mu    = zeros(1,factors);

go = 0;
while go==0
    try
        SIGMA = rand(factors);
        SIGMA(SIGMA==diag(SIGMA)) = 1;
        SIGMA = SIGMA - tril(SIGMA,-1) + triu(SIGMA,1)';
        r = mvnrnd(mu,SIGMA,N); go=1; % dirty way to test positive define matrix SIGMA
    end
end

parfor n=1:nboot
    data = mvnrnd(mu,SIGMA,N); 
    result        = rst_rep_anova_T2(data,gp,factors,0,[],[],0);
    robust_result = rst_rep_anova_T2(data,gp,factors,0,[],[],1); % default 20% trimmed mean
    P1(n) = result.p;
    P2(n) = robust_result.p;
end

type1_error =  mean([P1 P2]<0.05);


%% 3*3 sample repeated measure
% ----------------------------

N             = 200;                                        % 20 subjects
gp            = ones(N,1);                                 % 1 group
factors       = [3 2 3];                                   % repeated meaures

% test under H0, to check the type 1 error rate
nboot = 10000;
P1    = NaN(nboot,7);
P2    = NaN(nboot,7);

parfor n=1:nboot 
    data = randn(N,prod(factors)); 
    result        = rst_rep_anova_T2(data,gp,factors,0,[],[],0);
    robust_result = rst_rep_anova_T2(data,gp,factors,0,[],[],1); % default 20% trimmed mean
    P1(n,:) = result.p; 
    P2(n,:) = robust_result.p; 
end

type1_error =  [mean(P1<0.05)' mean(P2<0.05)'];

