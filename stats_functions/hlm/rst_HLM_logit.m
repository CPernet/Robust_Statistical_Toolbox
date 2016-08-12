% function results = rst_HLM_logit(data,N,regressor,groups,conditions)
% %
% % HLM for response curves using a hierarchical logistic linear model for 
% % binary data as for instance binary responses along a continuum or in time. 
% % Whilst Matlab default assumes common regression across all categories 
% % (parallel regression = proportional odds model) with different intercepts,
% % it is here always assumed that slope differ across conditions. 
% %
% % FORMAT: results = HLM_logit(data,N,regressor,groups,conditions)
% %
% % INPUT data       is a 3D matrix of dimension subjects, regressor, conditions
% %                  it represents the number of positive responses out of N
% %       N          is a scaler, vector or matrix for the total number of
% %                  events per subjects., sets and/or conditions if N is a scalar
% %                  the same N is used for all subjects, steps and
% %                  conditions, if N is a vector it represents the number of
% %                  events per subjects, if N is a matrix it must be 3D
% %                  matching data
% %       regressor  is a vector describing the sequence of events (eg time, steps, etc)
% %       groups     is a vector which describes the different groups of subjects
% %                  (e.g. [10 6] means 10 first subject are group 1 and 6 remaining 
% %                  subjects are group 2 ; set to [] if none)
% %       conditions is a vector which describes the different conditions
% %                  within subjects (e.g. 2 = 2 conditions but [2 2] = 2 factors 
% %                  with 2 levels each, ie 4 conditions)
% %
% % Analysis run 1st per subject.
% % 1. Response curves are vectorized and fitted to a design matrix
% % including a regressor (time or steps) per condition and a constant. 
% % This assumes data are ordinal, that is the regressor order has a meaning
% % and there are some commonality across condition to create a baseline.
% % 2. Fitted data are computed and coefficients are transformed for easier
% % interpretation at the result stage
% %
% % The group analysis
% % 1. Using the model parameters, a repeated measure ANOVA is set up testing
% % for the different effects - relies in Hotteling T^2
% % 2. Plots are produced using responses modelled at the 1st level.
% %
% % Requires the statistical toolbox
% % see mnrfit, mnrval
% % Cyril Pernet 16 OCtobre 2013
% 
% %% check inputs
% if nargin ~=4
%     error('wrong number of arguments in');
% end
% 
% % shift dimensions if needed
% if size(regressor,1) ~= 1;
%     regressor = regressor';
% end
% 
% if size(groups,1) ~= 1;
%     groups = groups';
% end
% 
% if size(conditions,1) ~= 1;
%     conditions = conditions';
% end
% 
% % check groups variable
% if isempty(groups)
%     groups = ones(size(data,1),1);
% end
% 
% % check N
%  if isscalar(N)
%      N = ones(size(data)).*N;
%  elseif isvector(N)
%      if size(N,1) == 1; N=N'; end
%      if size(N,1) ~= size(data,1)
%          error('the mumber of subjects in N doesn''t match the data');
%      end
%      N = repmat(N,[1,size(data,2),size(data,3)]);
%  end
%      
% % check correspondance between variables
% if numel(size(data))~=3
%     error('data must be a 3D matrix')
% elseif size(N) ~= size(data)
%     error('the mumber of events in the N variable doesn''t match the data');
% elseif length(regressor) ~= size(data,2)
%     error('the mumber of events in the regressor variable doesn''t match the data');
% elseif sum(groups) ~= size(data,1)
%     error('the mumber of subjects in the groups variable doesn''t match the data');
% elseif prod(conditions) ~= size(data,3)
%     error('the mumber of conditions in the conditions variable doesn''t match the data');
% end
% 
% 
% %% do the analysis per subject
% 
% % the design matrix
% X = [kron(eye(prod(conditions)),regressor')]; % ones are added by mnrfit
%  
%  for subject = 1:size(data,1)
%      % the data
%      Y1 = squeeze(data(subject,:,:));
%      Y2 = squeeze(N(subject,:,:));
%      prop = Y1./Y2;
%      % compute model
%      [B(subject,:),dev,stats]= mnrfit(X,[Y1(:) Y2(:)-Y1(:)],'model','ordinal','interactions','on','estdisp','on','link','logit');
%      % plot the fit
%      figure('Name',['Subject ' num2str(subject)]);
%      plot(regressor,prop,'o'); hold on
%      Yfit = mnrval(B(subject,:)',X);
%      Yhat{subject} = reshape(Yfit(:,1),size(Y1));
%      plot(regressor,Yhat{subject}); grid on
%      xlabel('regressor'); ylabel('odd ratio');
%      % get the slope at 0.5
%      for c=1:prod(conditions)
%          slope(subject,c) = B(subject,c+1)*0.05*[1-0.05];
%      end
%      title(['Slopes=' num2str(slope(subject,:))],'FontSize',12)
%  end
%  
%  
%  
%  %% do the group analysis
%  
%  % make groups a vector of labels
% if size(groups) == [size(data,1) 1]
%     GP = groups;
% else
%     GP = zeros(sum(groups),size(groups,2));
%     index1 = 1; index2 = groups(1);
%     for n=1:size(groups,2)
%         GP(index1:index2,n) = n;
%         if n<size(groups,2)
%             index1 = index1+groups(n);
%             index2 = index2+groups(n+1);
%         end
%     end
% end
% GP = sum(GP,2);
% 
% result = rst_rep_anova_T2(slope,GP,conditions,1000);
% 
%  %% return results and plots
%  
%  % compute the average response fit and robust 95%CI
%  %index1 = 1; index2 = groups(1);
%  %for g = 1:size(groups,2)
%   %   Y = Yhat{index1:index2}
%  
%   % plot the average residual fit
%  %stats{subject,winner}
% 
%  % plot the various effects from the ANOVA
%  %end
%  
%  
%  