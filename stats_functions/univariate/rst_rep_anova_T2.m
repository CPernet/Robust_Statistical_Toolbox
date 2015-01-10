function result = rst_rep_anova_T2(varargin)

% This function computes repeated measures ANOVAs - Unlike standard ANOVAs,  
% it relies on a multivariate framework which automatically account for the 
% correlation between measures. One advantage of this approach is that one  
% does not have to account for sphericity post hoc. In short we simply 
% either run a T2 on the repeated measures or a MANOVA (generalized T2 on 
% transformed data). The code implements equations described in Rencher 
% (2002) Methods of multivariate analysis John Wiley.
%
% FORMAT 
% result = rab_rep_anova_T2(data,gp,factors,nboot)
% result = rab_rep_anova_T2(data,gp,factors,nboot,{'name1','name2','etc'})
% result = rab_rep_anova_T2(data,gp,factors,nboot,{'name1','name2','etc'})
% result = rab_rep_anova_T2(data,gp,factors,nboot,{'name1','name2','etc'},C)
%
% INPUT
% Data:     2D matrix of subjects X factor level. Levels for all factors are
%           collapsed across columns in a hierarchical fashion. For example, a
%           a study with 2 factors, A and B, each with 2 levels (e.g., A1 and
%           A2), might organize data across columns like this: 
%           [A1B1   A1B2   A2B1   A2B2]
%           In this case, factor A is at the top of the hierarchy and factor B
%           is at the bottom. Therefore, factors A and B are the first and
%           second factors considered in the OUTPUT; e.g., result.F(1) for A
%           and result.F(2) for B
%
% gp:       Vector matched in length to rows of Data. Values indicate
%           level of a single between-group factor. Designs with multiple
%           between-group factors are not supported. For designs with no
%           between-group factors, set gp to []. 
%
% factors:	is a vector indicating the number of levels for each 
%           within-subject factor.
%           For example, [2 3] indicates 2 levels in factor A and 3 levels
%           in factor B. See info about Data to understand what is meant by
%           "factor A" and "factor B". In this case, there would have to be
%           2*3=6 columns in Data.
%
% nboot allows to run a bootstrap under H0 to estimate the distrubution of F 
% values - % WARNING if you don't have enough subjects, the resampling may 
% lead to singular or close to singular matrices, making the T2 invalid and 
% the boostrap F spurious. For no bootstap nboot=0
%
% 'name' default are A,B,C,AB,AC .. but you can input your owns if you wish 
% as a cell array
%
% C is optional to run contrasts on specified columns
%
%
% OUTPUT
%
% result.F F value of the Hotelling T2 test
% result.p corresponding level of significance
% result.names attribute a letter per factor by default
% result.pboot p value obtained using bootstrap (optional)
%
% EXAMPLE result = rep_anova_T2(data,gp,[2 2],0)
%         data   = [1 2 3 8 5 2 4 8 7 2 4; 1 5 7 9 4 2 6 4 8 3 5; ...
%                   4 6 8 4 1 8 0 1 4 6 2; 8 9 6 8 5 8 2 7 5 9 4]';  
%           gp   = [1 1 1 1 1 1 2 2 2 2 2]'; % 2 gps with 6 and 5 subjects
%
%      The data are organized as usual, like with SPSS, Statistica, etc ..
%
%      GP    factor A     level 1        level 2
%            factor B level 1 level 2 level 1 level 2 
%      1                1       1       4       8
%      1                2       5       6       9
%      1                3       7       8       6
%      1                8       9       4       8
%      1                5       4       1       5
%      1                2       2       8       8
%      2                4       6       0       2
%      2                8       4       1       7
%      2                7       8       4       5
%      2                2       3       6       9
%      2                4       5       2       4
%              
% Cyril Pernet V1 April 2011
% Carl M. Gaspar & GAR Dec 2011: edited the help comments.
% C. Pernet - updated January 2012 (thanks to Rik Henson for 
% making me adding cell names + pointing a possible error in df 
% for special cases including > 4 factors with 2 levels each 
% (still in debate)) 
% GAR June 2013: fixed bootstrap p value - bug spotted by Gil Shlomo
% Sharvit & Corrado Corradi Dell'Acqua

%% input stuff
% -------------
if nargin < 4 || nargin > 6
    error('wrong number of arguments')

else
    
    Data = varargin{1};
    gp   = varargin{2};
    if size(gp,1) ~= size(Data,1)
        error('data and gp parameter are of different size')
    elseif size(gp,2) > 1
        error('Designs with more than 1 between factor are not supported')
    end
    
    % deal with group structure
    % -------------------------
    if isempty(gp)
        gp = ones(size(Data,1), 1);
    else
        gp_values = unique(gp);
        k         = length(gp_values);
        
        % build the design matrix of gps
        X = NaN(size(gp,1),k+1);
        for g =1:k
            X(:,g) = gp == gp_values(g);
        end
        X(:,end) = 1;
    end
    
    % other parameters
    % ---------------
    factors = varargin{3};
    nboot = varargin{4};
    if ~isnumeric(nboot)
        error('nboot argument must be numeric')
    elseif nboot < 0
        nboot = abs(nboot); disp('nboot was negative, changed to positive')
    elseif nboot > 0
        Data_centered = NaN(size(Data));
        for g =1:k
            index = find(gp == gp_values(g));
            Data_centered(index,:) = detrend(Data(index,:),'constant');
        end
    end
    mynames = [];
    C = [];

    if nargin == 6
        C = varargin{6};
        mynames = varargin{5};
        if isnumeric(varargin{5})
            error('names can''t be numerics')
        end
    elseif nargin == 5
        mynames = varargin{5};
        if isnumeric(varargin{5})
            error('names can''t be numerics')
        end
    end
    
    if ~isempty(mynames)
        if size(mynames,2) ~= length(factors) || ~iscell(mynames)
            mynames = []; disp('names don''t match, switching to default');
        end
    end
end


%% basic info about the design
% -----------------------------
[n,p]         = size(Data);
nb_factors    = size(factors,2);
nb_effects    = (2^nb_factors - 1);
nb_conditions = prod(factors);

if unique(gp) == 1
    % one sample
    if nb_factors ==1
        type = 1;
    elseif nb_factors >1
        type = 2;
    end
else
    % k samples
    if nb_factors ==1
        type = 3;
    elseif nb_factors >1
        type = 4;
    end
end

%% analyze
% ---------
switch type
    
% One sample repeated measure
% ---------------------------
 
    % ---------------------------------------------------------------------
    case{1}  % 1 factor
        % -----------------------------------------------------------------
        if isempty(C)
            C = [eye(p-1) ones(p-1,1).*-1];  % contrast matrix
        end
        S = cov(Data); % covariance to account for spericity
        
        df           = p-1; 
        dfe          = n-p+1; 
        y            = nanmean(Data,1)'; % these are the means to compare
        Tsquare = n*(C*y)'*inv(C*S*C')*(C*y); % Hotelling Tsquare
        result.F     = ( dfe / ((n-1)*df) ) * Tsquare;
        result.p     = 1 - fcdf(result.F, df, dfe);
        result.df    = [df dfe];
        
        % bootstrap under H0
        % ------------------
        if nboot > 0
            dof_ratio = ( dfe / ((n-1)*df) );
            table = ceil(rand(size(Data,1),nboot).*size(Data,1));
            for B=1:nboot
                fprintf('bootstrap %g\n',B)
                boot_Data = Data_centered(table(:,B),:);
                S = cov(boot_Data);
                y = nanmean(boot_Data,1)';
                Tsquare = n*(C*y)'*inv(C*S*C')*(C*y);
                F(B) = dof_ratio * Tsquare;
            end
            % get p value
            v = sum(F < result.F)/nboot; % = how many times the observed F value > F under H0
            result.pboot =2*min(v,1-v); 
        end

        
        % -----------------------------------------------------------------
    case{2} % many factors
        % -----------------------------------------------------------------
        if isempty(C)
            [C,Names] = OrthogContrasts(factors); % set of orthogonal contrasts between factors
        end
        S = cov(Data); % covariance to account for spericity
        
        y = nanmean(Data,1)'; % these are the means to compare   
        if iscell(C)
            for effect = 1:size(C,2)
                c                = C{effect};
                df(effect)       = rank(c);
                dfe(effect)      = n-df(effect);
                Tsquare          = n*(c*y)'*inv(c*S*c')*(c*y); % is also t = sqrt(n*(c*y)'*inv(c*S*c')*(c*y));
                result.F(effect) = ( dfe(effect) / ((n-1)*(df(effect))) ) * Tsquare;
            end
        else % contrast
            Names    = [];
            df       = rank(C);
            dfe      = n-df;
            Tsquare  = n*(C*y)'*inv(C*S*C')*(C*y); 
            result.F = ( dfe / ((n-1)*(df)) ) * Tsquare;
        end
        result.p                      =  1 - fcdf(result.F, df, dfe);
        result.repeated_measure.names = Names;
        result.df                     = [df ; dfe];
        
        % bootstrap under H0
        % ------------------
        if nboot > 0
            table = ceil(rand(size(Data,1),nboot).*size(Data,1));
            for B=1:nboot
                fprintf('bootstrap %g\n',B)
                boot_Data = Data_centered(table(:,B),:);
                S = cov(boot_Data);
                y = nanmean(boot_Data,1)';
                if iscell(C)
                    for effect = 1:size(C,2)
                        c = C{effect};
                        Tsquare = n*(c*y)'*inv(c*S*c')*(c*y); % is also t = sqrt(n*(c*y)'*inv(c*S*c')*(c*y));
                        F(effect,B) = ( dfe(effect) / ((n-1)*(df(effect))) ) * Tsquare;
                    end
                else
                    Tsquare = n*(C*y)'*inv(C*S*C')*(C*y); 
                    F(B) = ( dfe / ((n-1)*(df)) ) * Tsquare;
                end
            end
            
            % get p values
            if iscell(C)
                for e = 1:nb_effects
                    v = sum(F(e,:) < result.F(e))/nboot;
                    result.pboot(e) = 2*min(v,1-v);
                end
            else
                v = sum(F < result.F)/nboot;
                result.pboot = 2*min(v,1-v);
            end
        end

        
% k samples repeated measure
% ---------------------------

        % -----------------------------------------------------------------
    case{3} % 1 within and 1 between factors
        % -----------------------------------------------------------------
        
        % get the error matrix of the full model
        % ---------------------------------------------------
        R  = eye(size(Data,1)) - (X*pinv(X));
        E  = (Data'*R*Data);                                                         
                       
        
        % compute the repeated measure
        % ------------------------------
        if isempty(C)
            C = [eye(p-1) ones(p-1,1).*-1]; % orthogonal contrast between levels
        end
        
        df  = p-1; 
        dfe = (n-p+1) - (k-1); % remove from dfe nb_gp - 1
        yp  = nanmean(Data,1)'; % average across gp
        ve  = sum(sum(X(:,1:end-1))-1); % - rank(X); % dfe for different sample sizes (gives the same as rank(X)*(sum(X(:,1))-1) for equal sample sizes            
        Spl = E/ve; % covariance of data split per gp 
        Tsquare  = n*(C*yp)'*inv(C*Spl*C')*(C*yp);
        result.repeated_measure.F    = ( dfe / (ve*df) ) * Tsquare;
        result.repeated_measure.p    = 1 - fcdf(result.repeated_measure.F, df, dfe);
        result.repeated_measure.df   = [df; dfe];
        
        % compute the gp effect ( = univariate between gps)
        % -------------------------------------------------            
        Y  = mean(Data,2); % average repeated measures
        [result.gp.F,result.gp.p,df,dfe] = local_glm(Y,X,k,sum(X(:,1:k),1),1);
        result.gp.df = [df,dfe];

       % compute the interaction (= multivariate on differences)
       % -------------------------------------------------------
       I        = (C*Data')';
       [result.interaction.F, result.interaction.p,df,dfe] = local_glm(I,X,k,sum(X(:,1:k),1),2);
       result.interaction.df = [df dfe];
       
        % bootstrap under H0
        % ------------------
        if nboot > 0
            table = ceil(rand(size(Data,1),nboot).*size(Data,1));
            for B=1:nboot
                fprintf('bootstrap %g\n',B)
                boot_Data = Data_centered(table(:,B),:);
                R  = eye(size(boot_Data,1)) - (X*pinv(X));
                E  = (boot_Data'*R*boot_Data);
                % compute the repeated measure
                yp  = mean(boot_Data,1)'; 
                Spl = E/ve;
                Tsquare = n*(C*yp)'*inv(C*Spl*C')*(C*yp);
                repeated_measure_F(B) = ( dfe / (ve*df) ) * Tsquare;
                % compute the gp effect
                Y = mean(boot_Data,2);
                [gp_F(B),p] = local_glm(Y,X,k,sum(X(:,1:k),1),1);
                % compute the interaction 
                I = (C*boot_Data')';
                [interaction_F(B), p] = local_glm(I,X,k,sum(X(:,1:k),1),2);
            end
            % get p values
            v = sum(repeated_measure_F < result.repeated_measure.F)/nboot;
            result.repeated_measure.pboot = 2*min(v,1-v);
            v = sum(gp_F < result.gp.F)/nboot;
            result.gp.pboot = 2*min(v,1-v);
            v = sum(interaction_F < result.interaction.F)/nboot;
            result.interaction.pboot = 2*min(v,1-v);
        end
        
        % -----------------------------------------------------------------
    case{4} % several within and 1 between factors
        % -----------------------------------------------------------------
                        
        % get the effect and error matrices of the full model
        % ---------------------------------------------------
        R  = eye(size(Data,1)) - (X*pinv(X));
        E  = (Data'*R*Data);                                                                              
        
        % compute the repeated measure
        % ------------------------------
        if isempty(C)
            [C,Names] = OrthogContrasts(factors); % set of orthogonal contrasts between factors
        end
        
        y  = mean(Data,1)'; % average across gp
        ve = 0; % dfe as a function of the numnber of subjects per gp
        for g=1:k
            v = sum(sum(X(:,g)==1));
            ve = ve + (v-1);
        end
        
        if iscell(C)
            for effect = 1:length(C)
                c           = C{effect};
                df(effect)  = rank(c);
                dfe(effect) = n-df(effect)-(k-1);
                Spl         = E/ve;
                Tsquare     = n*(c*y)'*inv(c*Spl*c')*(c*y);
                result.repeated_measure.F(effect) = ( dfe(effect) / (ve*df(effect)) ) * Tsquare;
                result.repeated_measure.p(effect) = 1 - fcdf(result.repeated_measure.F(effect), df(effect), dfe(effect));
            end
        else % contrast
            Names   = [];
            df      = rank(C);
            dfe     = n-df-(k-1);
            Spl     = E/ve;
            Tsquare = n*(C*y)'*inv(C*Spl*C')*(C*y);
            result.repeated_measure.F = ( dfe / (ve*df) ) * Tsquare;
            result.repeated_measure.p = 1 - fcdf(result.repeated_measure.F, df, dfe);
        end
        result.repeated_measure.names  = Names;
        result.repeated_measure.df     = [df ; dfe];
        
        % compute the gp effect (=univariate stat)
        % ---------------------------------------     
        Y  = mean(Data,2); % average repeated measures
        [result.gp.F, result.gp.p, df,dfe] = local_glm(Y,X,k,sum(X(:,1:k),1),1);
        result.gp.df = [df dfe];
        
       % compute the interactions with gp
       % ---------------------------------              
       if iscell(C)
           for effect = 1:length(C)
               c = C{effect};
               I = (c*Data')';
               [result.interaction.F(effect), result.interaction.p(effect),df,dfe]= local_glm(I,X,k,sum(X(:,1:k),1),2);
               Interaction_names{effect} = sprintf('Gp x %s',Names{effect});
               dof(:,effect) = [df dfe];
           end
       else
           I = (C*Data')';
           [result.interaction.F, result.interaction.p,df,dfe]= local_glm(I,X,k,sum(X(:,1:k),1),2);
           Interaction_names = 'Gp interaction';
           dof = [df dfe];
       end
       result.interaction.names = Interaction_names;
       result.interaction.df    = dof;
       
       % bootstrap under H0
       % ------------------
       if nboot > 0
           table = ceil(rand(size(Data,1),nboot).*size(Data,1));
           for B=1:nboot
               fprintf('bootstrap %g\n',B)
               boot_Data = Data_centered(table(:,B),:);
               R  = eye(size(Data,1)) - (X*pinv(X));
               E  = (Data'*R*Data);
               % compute the repeated measure
               y  = mean(boot_Data,1)'; 
               ve = 0; 
               for g=1:k
                   v = sum(sum(X(:,g)==1));
                   ve = ve + (v-1);
               end
               
               if iscell(C)
                   for effect = 1:length(C)
                       c   = C{effect};
                       df  = rank(c);
                       dfe = n-df-(k-1);
                       Spl = E/ve;
                       Tsquare = n*(c*y)'*inv(c*Spl*c')*(c*y);
                       repeated_measure_F(effect,B) = ( dfe / (ve*df) ) * Tsquare;
                   end
               else
                   df  = rank(C);
                   dfe = n-df-(k-1);
                   Spl = E/ve;
                   Tsquare = n*(C*y)'*inv(C*Spl*C')*(C*y);
                   repeated_measure_F(B) = ( dfe / (ve*df) ) * Tsquare;
               end
               
               % compute the gp effect (=univariate stat)
               Y = mean(boot_Data,2); 
               [gp_F(B), p] = local_glm(Y,X,k,sum(X(:,1:k),1),1);
               
               % compute the interactions with gp
               if iscell(C)
                   for effect = 1:length(C)
                       c = C{effect};
                       I = (c*boot_Data')';
                       [interaction_F(effect,B), p]= local_glm(I,X,k,sum(X(:,1:k),1),2);
                   end
               else
                   I = (C*boot_Data')';
                   [interaction_F(B), p]= local_glm(I,X,k,sum(X(:,1:k),1),2);
               end
           end
           
            % get p values
            if iscell(C)
%                 pval<-sum(test$TEST<=testb)/nboot
                for e = 1:length(C)
                    result.repeated_measure.pboot(e)=sum(repeated_measure_F(e,:)>=result.repeated_measure.F(e))./nboot;
                    result.interaction.pboot(e)=sum(interaction_F(e,:)>=result.interaction.F(e))./nboot;
                end
            else
                result.repeated_measure.pboot=sum(repeated_measure_F>=result.repeated_measure.F)/nboot;
                result.interaction.pboot=sum(interaction_F>=result.interaction.F)./nboot;
            end
            result.gp.pboot=sum(gp_F>=result.gp.F)./nboot;
       end
        
% --
end % closes the switch

% update names
if ~isempty(mynames)
for n=1:nb_factors; Names{n} = mynames{n}; end

index = nb_factors+1; for n=2:nb_factors
    interaction = nchoosek([1:nb_factors],n);
    for m=1:size(interaction,1)
        tmp = [Names{interaction(m,1)} ' by ' Names{interaction(m,2)}];
        if size(interaction,2) > 2
            for o = 3:size(interaction,2)
                tmp = [tmp ' by ' Names{interaction(m,o)}];
            end
        end
        Names{index} = tmp;
        index = index + 1;
    end
end

result.repeated_measure.names = Names;
if isfield(result,'interaction')
    result.interaction.names = Names;
end
end


end % closes the main function 

%% subfunctons

% ----------------------------------------------------
function [F,p,df,dfe] = local_glm(Y,X,nb_gp,nb_subjects,flag)

% compute the F and p values of a 1 way ANOVA - this is used here to 
% compute the effect of the gp in the repeated measure ANOVA or the 
% interactions - interactions are like testing a MANOVA on gps except
% the data are now differences between repeated measures
% the code works using a projection matrix, implementing the equations
% describe in Christensen, R. 2002. Plane answers to complex questions. 3rd
% Ed. Springer-Verlag
%
% FORMAT
% [F, p]= local_glm(Y,X,nb_gp,nb_subjects,flag)
%
% INPUTS
% Y the data
% X the design matrix
% nb_gp the number of independent groups
% nb_subjects the number of subjects per group
% flag = 2 for Hotelling stat otherwise univarite results are computed 
%
% OUTPUTS
% F the univariate or multivariate (Hotelling) F value
% p the corresponing level of significance


T        = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));  % SS Total
R        = eye(size(Y,1)) - (X*pinv(X));  % Residual matrix
E        = (Y'*R*Y);   % SS Error
Betas    = pinv(X)*Y;    
C = eye(size(X,2));
C(:,size(X,2)) = 0;
C0   = eye(size(X,2)) - C*pinv(C);
X0   = X*C0;                                                              
R0   = eye(size(Y,1)) - (X0*pinv(X0));
M    = R0 - R;    % M is the projection matrix onto Xc
H    = (Betas'*X'*M*X*Betas);    % SS Hypothesis (Effect)
df   = rank(X)-1;
dfe  = size(Y,1)-rank(X);

if flag == 2
    Eigen_values = decomp(E,H);
    p = size(Y,2); % = number of variables (dimension)
    vh = nb_gp - 1; % df = q above
    s = min(vh,p); % subspace in which mean Ys are located
    if sum(nb_subjects == nb_subjects(1)) == length(nb_subjects)
        ve = nb_gp*(nb_subjects(1)-1);     % dfe equal sample sizes
    else
        ve = sum(nb_subjects) - nb_gp;     % dfe different sample sizes
    end
    
    if s > 1
        m = (abs(vh-p)-1)/2;
        N = (ve-p-1) / 2;
        U = sum(Eigen_values);
        df_Hotelling = s*(2*m+s+1); 
        dfe_Hotelling = 2*(s*N+1);
        F = (dfe_Hotelling*U) / (s^2*(2*m+s+1));
        p = 1-fcdf(F, df_Hotelling, dfe_Hotelling);
        df = df_Hotelling; dfe = dfe_Hotelling;
    else % = only one non zeros Eigen value s = 1 and/or vh = 1
        U = max(Eigen_values);
        df_Hotelling = p;
        dfe_Hotelling = ve-p+1;
        F = (dfe_Hotelling/df_Hotelling) * max(Eigen_values);
        p = 1-fcdf(F,df_Hotelling,dfe_Hotelling);
    end
else
    F    = (diag(H)./(rank(X)-1)) ./ (diag(E)/(size(Y,1)-rank(X)));
    p    = 1 - fcdf(F, (rank(X)-1), (size(Y,1)-rank(X)));
end
end

% ------------------------------------------
function [M,MList]= OrthogContrasts(nTreats)

% Generates coefficients corresponding to one particular set of
% orthogonal contrasts for a set of variables
%
% FORMAT
% [contrasts, names]= OrthogContrasts([2 2])
%
% INPUT
% nTreats is a vector with the levels of the various factors in an ANOVA
%
% OUTPUT
% M is a cell with the various contrast in it
% Mlist is a series of set names
%
% Taken from Matthew Nelson GenOrthogComps.m
% simplified by Cyril Pernet 16-09-2009

if length(nTreats) == 1
    error('there must be at least 2 factors')
end

alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ';

dfs=nTreats-1;
nFacs=length(nTreats);
nCombs=prod(nTreats);

% Conpute main effects contrasts
M=cell(1,nFacs);
for iFac=1:nFacs
    M{iFac}=repmat(0,dfs(iFac),nCombs);
    
    if iFac==nFacs;
        nLowCombs=1;
    else
        nLowCombs=prod(nTreats(iFac+1:end));
    end
    
    if iFac==1
        nUpCombs=1;
    else
        nUpCombs=prod(nTreats(1:iFac-1));
    end
    
    for idf=1:dfs(iFac)
        tmpComps=repmat( 0,1,nLowCombs*nTreats(iFac) );
        tmpComps( 1:nLowCombs*idf )=1/idf;
        tmpComps( nLowCombs*idf+1:nLowCombs*(idf+1) )=-1;
        tmpComps=repmat(tmpComps,1,nUpCombs);
        M{iFac}(idf, : )=tmpComps;
        MList{iFac}=alphabet(iFac);
    end       
end

% Do interactions as products of the main effect contrasts
if nFacs>1        
    for curnFacs=2:nFacs     %start with 2 way interactions and go up to nFacs-way interactions
        for initFac=1:nFacs-(curnFacs-1)                       
            [M MList]=RLoop(curnFacs,initFac,repmat(1,1,nCombs), length(M)+1, '', M,MList,dfs,alphabet); 
        end
    end
end
end

% ------------------------------------------------------------------------------------
function [M MList]=RLoop(FacRem,curFac,curRow,curMNum, nextList, M,MList,dfs,alphabet) 

%given the initial inputs, and the initial curFac, this should calc all the
%downstream M vals for that initial curFac... routine called by OrthogContrasts

FacRem=FacRem-1;
for idf=1:dfs(curFac)
    if FacRem==0
        if curMNum>length(M);       
            M{curMNum}=[];      
            MList{curMNum}=[nextList alphabet(curFac)];
        end
        M{curMNum}(end+1,:)=curRow.*M{curFac}(idf,:);
    else
        nextRow=curRow.*M{curFac}(idf,:);
        nextList(end+1)=alphabet(curFac);
        baseFacNum=curFac+1;     %this is needed for the numbering of the outputs...
        for iNextFac=curFac+1:length(dfs)-FacRem+1
            [M MList]=RLoop( FacRem,iNextFac,nextRow,curMNum+(iNextFac-baseFacNum), nextList, M,MList,dfs,alphabet );
        end
    end
end
end

% ---------------------------------
function eigen_values = decomp(E,H)

% this function is used to decompose inv(E)*H
% Following Rencher 2002 (Methods of multivariate
% analysis - Wiley) we note that eig(inv(E)*H) =
% eig((E^1/2)*H*inv(E^1/2)) = eig(inv(U')*H*inv(U))
%
% E^1/2 is the square root matrix of E and U'U = E
% (Cholesky factorization). Using the Cholesky 
% factorisation, we return positve eigen values
% from inv(U')*H*inv(U) which is posite semidefinite.
% If this procedre fails (E is not positive definite) 
% we then use an SVD decomposition (the SVD decomposition
% is taken from the SPM code in spm_voi.m by Karl Friston
% Wellcome Trust Centre for Neuroimaging).


try
    U = chol(E);
    eigen_values = eig(inv(U')*H*inv(U));

catch
    y = (pinv(E)*H);
    [m n]   = size(y);
    if m > n
        [v s v] = svd(y*y');
        s       = diag(s);
        v       = v(:,1);
        u       = y*v/sqrt(s(1));
    else
        [u s u] = svd(y'*y);
        s       = diag(s);
        u       = u(:,1);
        v       = y'*u/sqrt(s(1));
    end
    d  = sign(sum(v)); u = u*d;
    eigen_values  = u*sqrt(s(1)/n);
end
end
