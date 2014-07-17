function varargout = rst_rep_anova_plot(Data,factors,gp,alpha)

% plots effects for a repeated measure anova
% plots include confidence intervals for within-subjects
%
% INPUT
% Data    = 2D matrix (n*p) of repeated measures
% gp      = vector (n*1) indicating gp belonging  - doesn't do nested gp
% factors = vector indicating the levels of each factors like [2 3 4] --> (prod(factors=p))
% alpha   = 5% (default) or 10%
%
% OUTPUT
% y       = means of the difference conditions
% CI      = confidence intervals of the means
% optional (see boot_multicompare.m)
% diff    = all pairwise differences between conditions
% CI_boot = percentile bootstrap confidence intervals of the differences
% p       = percentile bootstrap p values
% h       = significance of the p value (adjusted for multiple tests)

% CI are computed as described in Loftus, G.R (2002) Analysis,
% Interpretation and Visual Presentaion of Experimental Data.
% In Steven's handbook of experimental psychology (3rd Ed): 
% Methodology in experimental psychology (volume 4). 
% John Wiley & sons, Inc., NewYork
%
% Cyril Pernet 16-08-2011

%% basic info about the design
% -----------------------------
if nargin <4
    alpha = 5/100;
elseif nargin <3
    error('wrong number of arguments');
elseif nargin > 4
    error('wrong number of arguments');
end

[n,p]         = size(Data);
nb_factors    = size(factors,2);
nb_effects    = (2^nb_factors - 1);
nb_conditions = prod(factors);
nb_gp         = size(gp,2);
if nb_gp > 1
    errordlg('Designs with more than 1 between factor are not supported')
    return
end

%% analyze
% ---------
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

%% 
switch type
    
% One sample repeated measure
% ---------------------------
 
    % ---------------------------------------------------------------------
    case{1}  % 1 factor 
        % -----------------------------------------------------------------
        
        % CI based on the subject * factor interaction
        % --------------------------------------------
        
        % sum square total
        SST = (Data(:)-(repmat(mean(Data(:)),n*p,1)))'*(Data(:)-(repmat(mean(Data(:)),n*p,1)));
        
        % sum square factor
        C = [eye(p-1) ones(p-1,1).*-1];  % contrast matrix        
        % --------------------
        y = nanmean(Data,1)';
        varargout{1} = y;
        % ---------------------
        SSF = n*(C*y)'*(C*y);

        % sum square subject
        S = sum(Data,2);
        Y = Data(:);
        X = [kron(eye(p),ones(n,1)) repmat(S,p,1)];
        Betas    = pinv(X)*Y;
        R        = eye(size(Y,1)) - (X*pinv(X));
        C = zeros(size(X,2));
        C(end,end) = 1;
        C0   = eye(size(X,2)) - C*pinv(C);
        X0   = X*C0;                                                         
        R0   = eye(size(Y,1)) - (X0*pinv(X0));
        M    = R0 - R;                                                          
        SSS    = (Betas'*X'*M*X*Betas);                                         
        
        % sum square inteaction
        SSI = SST - SSF - SSS;
        dfI = length(Y) - (p-1) - (n-1) - 1;

        % plot
        t_critical= tinv(1-alpha/2,dfI);
        % --------------------
        CI = sqrt((SSI/dfI) / n) * t_critical; 
        varargout{2} = CI;
        % --------------------
        figure('Name', 'Repeated measure ANOVA')
        set(gcf,'Color','w')
        errorbar(y,repmat(CI,1,p),'r','LineWidth',2); 
        title('Factor 1','FontSize',16); grid on
        ylabel('data values','FontSize',14); 
        xlabel('conditions','FontSize',14);
        set(gca,'FontSize',12);
        set(gca,'xtick',[1:p])
        
        if nargout == 6
            [varargout{3},varargout{4},varargout{5},varargout{6}]= boot_multicompare(Data,gp,alpha);
        end

        % -----------------------------------------------------------------
    case{2} % many factors  
        % -----------------------------------------------------------------
        
        % CI based on pooled errors
        % ---------------------------
        
        % sum square total
        SST = (Data(:)-(repmat(mean(Data(:)),n*p,1)))'*(Data(:)-(repmat(mean(Data(:)),n*p,1)));
       
        S = sum(Data,2);
        Y = Data(:);
        X = [kron(eye(p),ones(n,1)) repmat(S,p,1)];
        Betas    = pinv(X)*Y;
        R        = eye(size(Y,1)) - (X*pinv(X));
        % all conditions
        C = eye(size(X,2));
        C(end,end) = 0;
        C0   = eye(size(X,2)) - C*pinv(C);
        X0   = X*C0;                                                         
        R0   = eye(size(Y,1)) - (X0*pinv(X0));
        M    = R0 - R;                                                          
        SSF    = (Betas'*X'*M*X*Betas);                                         
        % subjects
        C = zeros(size(X,2));
        C(end,end) = 1;
        C0   = eye(size(X,2)) - C*pinv(C);
        X0   = X*C0;                                                         
        R0   = eye(size(Y,1)) - (X0*pinv(X0));
        M    = R0 - R;                                                          
        SSS  = (Betas'*X'*M*X*Betas);                                         

        % pooled error
        SSI = SST - SSF - SSS;
        dfI = (nb_conditions-1)*(n-1);

        % plot
        t_critical= tinv(1-alpha/2,dfI);
        % --------------------
        CI = sqrt((SSI/dfI) / n) * t_critical;
        % --------------------
        q = eval(cell2mat(inputdlg(['which effect to plot? 1-->' num2str(nb_effects)],'choose effect'))); 
        tmp = nanmean(Data,1); 
        if q < nb_factors % main effects
            unit = p/prod(factors(1:q));
            tmp = mean(reshape(tmp,unit,length(tmp)/unit),1);
            unit = [1:factors(q)];
            increment = length(unit);
            y = zeros(1,factors(q));
            index = unit; roundx = 1;
            while index(end) <= length(tmp)
                if roundx==1
                    y = tmp(index);
                else
                    y = (y+tmp(index))/2;
                end
                index = index+increment;
                roundx = roundx+1;
            end
        elseif q == nb_factors 
            newdata = (reshape(tmp,factors(q),p/factors(q)))';
            y = mean(newdata);
        else % otherwise plot all
            y = tmp;
        end
        
        figure('Name', 'Repeated measure ANOVA')
        set(gcf,'Color','w')
        errorbar(y,repmat(CI,1,length(y)),'r','LineWidth',2);
        if q <= nb_factors % main effects
            title(['Factor ' num2str(q)],'FontSize',16); grid on
        else
            title('all conditions','FontSize',16); grid on
        end
        ylabel('data values','FontSize',14);
        xlabel('conditions','FontSize',14);
        set(gca,'FontSize',12);
        set(gca,'xtick',[1:p])

        % do multiple comparisons if 'possible'
        % ------------------------------------
        if nargout == 6
            multicomp = 1;
            if q < nb_factors % main effects
                unit = p/prod(factors(1:q));
                tmp = mean(reshape(tmp,unit,length(tmp)/unit),1);
                unit = [1:factors(q)];
                increment = length(unit);
                newData = zeros(1,factors(q));
                index = unit; roundx = 1;
                while index(end) <= length(tmp)
                    if roundx==1
                        newData = tmp(index);
                    else
                        newData = (y+tmp(index))/2;
                    end
                    index = index+increment;
                    roundx = roundx+1;
                end
            elseif q == nb_factors
                newdata = (reshape(tmp,factors(q),p/factors(q)))';
                newData = mean(newdata);
            elseif q ==  nb_effects
                newData = Data;
            else
                warndlg('multiple comparisons not performed - please average the relevant columns by hand and send these data to boot_multicompare.m [diff,CI_boot,p,h]= boot_multicompare(Data,gp,alpha)')
                multicomp = 0;
            end
            
            if multicomp ==1
                [varargout{3},varargout{4},varargout{5},varargout{6}]= boot_multicompare(newData,gp,alpha);
            end
        end

        
% k samples repeated measure
% ---------------------------

        % -----------------------------------------------------------------
    case{3} % 1 within and 1 between factors
        % -----------------------------------------------------------------
        
       
       
        % -----------------------------------------------------------------
    case{4} % several within and 1 between factors
        % -----------------------------------------------------------------
                        
% --
end % closes the switch
end % closes the main function 



