function varargout = rst_single_subject(rawAcc, rawRT, varargin)

% rab function to clean up and extract relevant information per subject
% 1st the raw data are cleaned using the Median Absolute Deviation (MAD)
% - a cutoff can be specified for RT and each condition is cleaned separately
% 2nd the mean, trimmed mean and median is extracted 
% optional: - figure output
%           - clean-up
%           - percentage of trimming 
%           - d prime
%           - vincentized RT        
%           - dynamic d prime
%
% [Acc,RT,EstAcc,EstRT] = bat_single_subject(rawAcc, rawRT);
% Acc = matrix of Accuracy data cleaned up
% RT = matrix of RT data cleaned up
% EstAcc = mean, trimmed mean and median accuracy
% EstRT = mean, trimmed mean and median RT
%
% [Acc,RT,EstAcc,EstRT,dprime,vincentRT,dynamicdprime] = ...
%        rab_single_subject(rawAcc, rawRT, 'figure','on', 'cleanup','on',[], ...
%        'percent',20, 'dprime','on','vintentizedRT','on','dynamicdprime','on');
% cleanup 'on' can be followed by [] or cutoff values e.g. [0.08 0.6]
% dprime = matrix of distance d' between all possible conditions
% vincentRT = cleaned RT sorted bins
% dynamicdprime = dprime matrices sorted on dim 3 per RT bins

% Cyril Pernet / Guillaume Rousselet 
% ----------------------------------
%  Copyright (C) RST Team 2014


% -------------
%% Input check
% -------------
if nargin < 2
    error('not enough arguments')
elseif nargin > 2
    if length(varargin) > 8
        error('too many options specified')
    else
        % defaults
        options.figure  = 'off';
        options.cleanup  = 'off';
        options.percent = 20/100;
        options.dprime = 'off';
        options.vintentizedRT = 'off';
        options.dynamicdprime = 'off';
        % check
        match = 0;
        for l=1:length(varargin);
            if strcmp('figure',varargin{l})
                options.figure = varargin{l+1}; match = 1;
            elseif strcmp('clean-up',varargin{l})
                options.cleanup = varargin{l+1}; match = 1;
                options.cleanupv = varargin{l+2}; 
            elseif strcmp('percent',varargin{l})
                options.percent = varargin{l+1}; match = 1;
            elseif strcmp('dprime',varargin{l})
                options.dprime = varargin{l+1}; match = 1;
            elseif strcmp('vintentizedRT',varargin{l})
                options.vintentizedRT = varargin{l+1}; match = 1;
            elseif strcmp('dynamicdprime',varargin{l})
                options.dynamicdprime = varargin{l+1}; match = 1;
            end
        end
        if match == 0
            error('option argument not recognized')
        end
    end
else
    options.percent = 20/100;
end

% check data dim
[n,p] = size(rawAcc);
if size(rawAcc) ~= size(rawRT)
    error('Accuracy and RT matrices have different dimensions')
end
quant = 0;
options

% -------------------------
%% Data clean up using MAD
% -------------------------
if isfield(options,'cleanup') == 1
    
    
    % 1st plot data 
    if isfield(options,'figure') == 'on'
        figure('Name','RT histograms')
        subplot(1,2,1); rab_density_hist(rawRT(:));
    end
    
    % 2nd overall cutoff
    % if isfield(options,'cleanupv')
    I = rab_madmedianrule(rawRT(:));
    lowbound = min(rawRT(find(I==0)));
    highbound = max(rawRT(find(I==0)));
    cutoff = cell2mat(inputdlg(['cutoff for RTs (' num2str(lowbound) ' ' num2str(highbound) ')'],'RTs cutoff'));
    if sum(cutoff) == 0
        if ischar(cutoff) % cutoff ='' if just pressed ok
            fprintf('using default values %g %g as low and high bounds \n',lowbound,highbound)
        else % cutoff = [] pressed cancel
            lowbound = min(rawRT(:));
            highbound = max(rawRT(:));
            fprintf('no cutoff applied')
        end
    else
        cutoff = str2num(cutoff);
        if cutoff(1) > cutoff(2)
            error('low bound bigger than high bound')
        else
            lowbound = cutoff(1);
            highbound = cutoff(2);
        end
    end
    good_trials = intersect(find(rawRT(:)>lowbound),find(rawRT(:)<highbound));
    subplot(1,2,2); rab_density_hist(rawRT(good_trials));
    
    % reset matrices with 1st cleanup
    tmpRT = NaN(size(rawRT));
    tmpRT(good_trials) = rawRT(good_trials);
    tmpAcc = NaN(size(rawAcc));
    tmpAcc(good_trials) = rawAcc(good_trials);
    
    % 3rd do a clean up per condition
    I = bat_madmedianrule(tmpRT);
    good_trials = find(I==0);
    cleanRT = NaN(size(rawRT));
    cleanRT(good_trials) = tmpRT(good_trials);
    varargout{1}.RT = cleanRT;
    cleanAcc = NaN(size(rawAcc));
    cleanAcc(good_trials) = tmpAcc(good_trials);
    varargout{2}.Accuracy = cleanAcc;
    
else
    varargout{1}.RT = RT; cleanRT = RT;
    varargout{2}.Accuracy = Acc; cleanAcc = Acc;
end

% ----------------
%% Get estimators 
% ----------------
percent = options.percent;
if percent > 1
    percent = percent / 100;
end

meanRT = nanmean(cleanRT,1);
trimmeanRT = nantrimmean(cleanRT,percent);
medianRT = nanmedian(cleanRT,1);
varargout{3}.summary_statRT.mean = meanRT;
varargout{3}.summary_statRT.trimmean = trimmeanRT;
varargout{3}.summary_statRT.median = medianRT;

meanAcc = nanmean(cleanAcc,1).*100;
trimmeanAcc = nantrimmean(cleanAcc,percent).*100;
medianAcc = nanmedian(cleanAcc,1).*100;
varargout{4}.summary_statAcc.mean = meanAcc;
varargout{4}.summary_statAcc.trimmean = trimmeanAcc;
varargout{4}.summary_statAcc.median = medianAcc;

% ---------
%% Options 
% ---------
index = 4;

% dprime - all pairs
% -------------------
if isfield(options,'dprime') == 1
    if strcmp('on',options.dprime) == 1
        combinations = nchoosek([1:p],2);
        dprime = zeros(p,p);
        for i=1:length(combinations)
            c1 = combinations(i,1);
            c2 = combinations(i,2);
            n1 = numel(find(~isnan(cleanAcc(:,c1)))); % the number of trial for stimulus 1
            n2 = numel(find(~isnan(cleanAcc(:,c2)))); % the number of trial for stimulus 2
            h = nansum(cleanAcc(:,c1)); % the number of Hits (yes on stimulus 1)
            f = nansum(cleanAcc(:,c2)); % the mumber of false alarms (yes on stimulus 2)
            [zHits,zFA,d,c,cprime,beta] = rab_sdt1(n1,n2,h,f,0);
            dprime(c1,c2) = d; dprime(c2,c1) = d;
        end
        index = index+1; varargout{index} = dprime;
    end
end

% Vincentized RT
% ---------------
if isfield(options,'vintentizedRT') == 1
    if strcmp('on',options.vintentizedRT) == 1
        quant = cell2mat(inputdlg('how many bins ?','Vincentization'));
        if sum(quant) == 0
            if ischar(quant) % cutoff ='' if just pressed ok
                disp('using default values of 3');
                quant = 3;
            else % pressed cancel
                disp('Vincentization cancelled')
            end
        else
            quant = str2num(quant);
        end
         
        if sum(quant) ~= 0
            for i=1:p
                tmp = cleanRT(:,i); tmp = tmp(find(~isnan(tmp)));
                distribution{i} = bat_vincentize(tmp,quant,0);
            end
        end
        index = index+1; varargout{index} = distribution;
    end
end

% Dynamic d'
% -----------
if isfield(options,'dynamicdprime') == 1
    if strcmp('on',options.dynamicdprime) == 1
        if quant ==0 % already used in Vincentization
            quant = cell2mat(inputdlg('how many bins ?','Vincentization'));
            if sum(quant) == 0
                if ischar(quant) % cutoff ='' if just pressed ok
                    disp('using default values of 3');
                    quant = 3;
                else % pressed cancel
                    disp('Vincentization cancelled')
                end
            else
                quant = str2num(quant);
            end
        end
        
        if sum(quant) ~= 0
            combinations = nchoosek([1:p],2);
            dynamic_dprime = zeros(p,p,quant);
            for i=1:length(combinations)
                c1 = combinations(i,1);
                tmp1 = cleanAcc(:,c1); tmp1 = tmp1(find(~isnan(tmp1)));
                tmp11 = cleanRT(:,c1); tmp11 = tmp11(find(~isnan(tmp11)));
                tmp111 = ones(length(tmp1),1);
                c2 = combinations(i,2);
                tmp2 = cleanAcc(:,c2); tmp2 = tmp2(find(~isnan(tmp2)));
                tmp22 = cleanRT(:,c2); tmp22 = tmp22(find(~isnan(tmp22)));
                tmp222 = ones(length(tmp2),1).*2;
                data = [[tmp11;tmp22] [tmp1;tmp2] [tmp111;tmp222]];
                out = bat_dynamicdprime(data,quant);
                for q=1:quant
                    dynamic_dprime(c1,c2,q) = out(q);
                    dynamic_dprime(c2,c1,q) = out(q);
                end
            end
        end
        index = index+1; varargout{index} = dynamic_dprime;
    end
end

pause; close
end % ends the function

% subfunction to compute the trimmed mean
function trimmean = nantrimmean(data,percent)
n = sum(~isnan(data),1);
trimming = round(percent.*n);
trimmean = NaN(1,size(data,2));
for i=1:size(data,2)
    tmp = sort(data(:,i)); tmp = tmp(find(~isnan(tmp)));
    low = trimming(i); high = n(i)-trimming(i);
    trimmean(i) = mean(tmp(low:high));
end
end
