function [t,p,CI]= two_ttest(gp1,gp2,type,aphav,graphout)

% 2 samples t-tests between grp1 and gp2 for different data structure type
% FORMAT [t,p,CI]= 2ttest(gp1,gp2,type)
% INPUT: gp1 and gp2 are two column vectors or matrices of data
%                    data can be of different sizes, use NaN for matrices
%        type can be 'equal' meaning equal variance is assumed
%                    'unequal' meaning unequal variance is assumned
%        alphav the the alpha value used to compute the confidence
%        intervals (default is 95%)
%        graphout 1/0 indicates graphical outputs (default is 0)
% OUTPUT t is the t value
%        p is the p-value
%        CI is the alphav confidence interval of the difference
%
% Cyril Pernet - 03 Jan 2015

%% input checks
if nargin == 3
    alphav = 5/100;
    graphout = 0;
elseif nargin == 4
    graphout = 0;
end

if nargout == 0
    graphout = 1;
end

if size(gp1,2) ~= size(gp2,2)
    error('different number of variables were found in gp1 and gp2, revise inputs')
end

%% compute
Difference = nanmean(gp1,1) - nanmean(gp2,1);

switch type % http://en.wikipedia.org/wiki/Student%27s_t-test#Independent_.28unpaired.29_samples

case('equal') % Student's t=test
    
    % if no NaN compute for matrices
    if sum(isnan(gp1(:)))+sum(isnan(gp2(:))) == 0
        if size(gp1,1) == size(gp2,1) % equal sample sizes
            s12 = sqrt((var(gp1,0,1)+var(gp2,0,1))./2);
            se = (s12*sqrt(2/size(gp1,1)));
            df = 2*size(gp1,1) - 2 ;
        else % unequal sample sizes
            A = (size(gp1,1)-1)*nanvar(gp1,0,1)+(size(gp2,1)-1)*nanvar(gp2,0,1);
            df = size(gp1,1)+size(gp2,1)-2; s12 = sqrt(A/df); 
            se = (s12*sqrt((1/size(gp1,1))+(1/size(gp2,1))));
        end
        t = Difference / se;
    else
        parfor column = 1:size(gp1,2)
            n1=sum(~isnan(gp1(:,column)));
            n2=sum(~isnan(gp2(:,column)));
            if n1 == n2
                s12 = sqrt((nanvar(gp1(:,column),0,1)+nanvar(gp2(:,column),0,1))./2);
                df = 2*n1 - 2 ; se(column) = s12*sqrt(2/n1);
            else % unequal sample sizes, equal variance
                A = (n1-1)*nanvar(gp1(:,column),0,1)+(n2-1)*nanvar(gp2(:,column),0,1);
                df = n1+n2-2; s12 = sqrt(A/df); se(column) = s12*sqrt((1/n1)+(1/n2));
            end
            t(column) = Difference(column) / se(column);
        end
    end
    
case('unequal') % Welch's t-test
    
    if sum(isnan(gp1(:)))+sum(isnan(gp2(:))) == 0
        n1 = size(gp1,1); n2 = size(gp2,1);
        A = (var(gp1,0,1)./n1)+(var(gp2,0,1)./n2);
        se = sqrt(A); t = Difference / se;
        B = ((var(gp1,0,1)/n1)^2)/(n1-1);
        C = ((var(gp2,0,1)/n2)^2)/(n2-1);
        df = A^2 / (B+C); % Welch–Satterthwaite 
    else
        parfor column=1:size(gp1,2)
            n1=sum(~isnan(gp1(:,column)));
            n2=sum(~isnan(gp2(:,column)));
            A = (nanvar(gp1(:,column),0,1)/n1)+(nanvar(gp2(:,column),0,1)/n2);
            B = ((nanvar(gp1(:,column),0,1)/n1)^2)/(n1-1);
            C = ((nanvar(gp2(:,column),0,1)/n2)^2)/(n2-1);
            se(column) = sqrt(A); df = A^2 / (B+C); % Welch–Satterthwaite 
            t(column) = Difference / se(column);
        end
    end

otherwise
   disp('Unknown type')
end

p = 2 * tcdf(-abs(t),df);
if nargout > 2
    spread = tinv(1 - alphav ./ 2, df) .* se;
    CI = [Difference-spread, Difference+spread];
end

%% graphical output
if graphout == 1
    
end