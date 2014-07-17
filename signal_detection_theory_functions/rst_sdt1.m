function [zHits,zFA,dprime,c,cprime,beta] = rst_sdt1(varargin)
% implement signal detection theory
% based on the book Detection theory A user guide
% Authors N.A. Macmillan and C.D. Creelman
% one-interval design - two-stimuli needs to be discriminated
%
% [zHits,zFA,dprime,c,cprime,beta] = bat_sdt1(100,100,80,30,1);
%
% INPUT
% 1st value is the number of trial for stimulus 1 
% 2st value is the number of trial for stimulus 2 
% 3nd value is the number of Hits (yes on stimulus 1)
% 4rd value is the mumber of false alarms (yes on stimulus 2)
% 5th value indicates to display results (any value) or not (0)
%
% OUTPUT
% zHits, zFA, dprime = prob. of hit and FA and there diff.
% c and cprime = criterion location and relative criterion location
% beta = likelihhod ratio
%
% Cyril Pernet Jan 2011.

if isempty(varargin)
    
    disp(' ');
    disp('--------------------------------------------------');
    disp('Stimulus class               Response             ');
    disp('              ------------------------------------');
    disp('                   YES                NO          ');
    disp('              ------------------------------------');
    disp('Stimulus 1         Hits             Misses        ');
    disp('Stimulus 2     False alarms   Correct rejection   ');
    disp('--------------------------------------------------');
    disp(' ');
    
    disp('The number of stimululi is assumed equal for S1 and S2');
    N1       = input('how many trials for stimulus1?   ');
    N2       = input('how many trials for stimulus2?   ');
    Hits     = input('how many Hits?                 ');
    FA       = input('how many False alarms?         ');
    makeplot = 1;
    
else
    N1       = varargin{1};
    N2       = varargin{2};
    Hits     = varargin{3};
    FA       = varargin{4};
    makeplot = varargin{5};
end

% basic stuff so that it fits the description
Misses = N1 - Hits;
CR     = N2 - FA;

if makeplot ~= 0
    disp   (' ');
    disp   ('--------------------------------------------------');
    disp   ('Stimulus class               Response             ');
    disp   ('              ------------------------------------');
    disp   ('                   YES                NO          ');
    disp   ('              ------------------------------------');
    fprintf('Stimulus 1         %g                 %g          ',Hits,Misses);
    disp   (' ');
    fprintf('Stimulus 2         %g                 %g          ',FA,CR);
    disp   (' ');
    disp   ('--------------------------------------------------');
    fprintf('%g trials stimulus 1, %g trials stimulus 2',N1,N2);
    disp   (' ');
end

% compute P(H) and P(F)
% if we get 1 or 0 we choose a value 1/(2*N);

if (Hits/N1)  == 1
    Hits = 1-(1/(2*N1));
elseif Hits == 0
    Hits = 1/(2*N1);
else
    Hits = Hits/N1;
end

if (FA/N2)  == 1
    FA = 1-(1/(2*N2));
elseif FA == 0
    FA = 1/(2*N2);
else
    FA = FA/N2;
end

% compute the z values of P(H) and P(F) and d'

zHits = norminv(Hits,0,1);
zFA   = norminv(FA,0,1);
dprime = (zHits) - (zFA);

% compute the criteria and response bias
% note that a positive bias is a tendency to say 'no'
% and a negative bias a tendency to say 'yes'

% criterion location
temp = (zHits) + (zFA);
c    = (-1/2) .* (temp);

% relative c
cprime = (c)./(dprime);

% likelihhod ratio beta
beta = exp(c.*dprime);

% display results
if makeplot ~= 0
    
    for i = 1:length(zHits)
        disp(' ');
        fprintf('z(Hits)= %g ; z(False alarms) = %g ; d''= %g',zHits(i),zFA(i),dprime(i));
        disp(' ');
        fprintf('criterion location= %g ; relative c= %g ; likelihood ratio= %g',c(i),cprime(i),beta(i));
        disp(' ');
    end
    
    zROC(zFA,zHits)
end
end

% embeded function to make the zROC plot
function zROC(zFA,zHits)

x = [-3:1:3];
y = [-3:1:3];
plot(x,y,'r');
hold on
plot(zFA,zHits,'*');
hold off
grid
title('zROC');
xlabel('z P(False alarms)');
ylabel('z P(Hits)');
end
