function [tx,g]=rst_tsample(x,percent)

% function [tx]=tsample(x,percent)
% Returns the trimmed sample of x.
% x must be a vector
% percent must be a number between 0 and 100
% The trimmed sample is calculated by dropping the lower and upper g values from
% x, where g=floor((percent/100)*length(x)). If x is empty, then NaN is
% returned.
%
% Original code provided by Prof. Patrick J. Bennett, McMaster University
%
% See also TM

% The output size for [] is a special case, handle it here.
if isequal(x,[]), tx = NaN; return; end;

% if nargin < 2
%     error('tsample:TooFewInputs', 'tsample requires two input arguments.');
% elseif percent >= 100 || percent < 0
%     error('tsample:InvalidPercent', 'PERCENT must be between 0 and 100.');
% end

% make sure that x is a vector
sz = size(x);
if sz > 2 | min(sz) > 1    
  error('tsample requires x to be a vector, not a matrix.');
end

n = length(x);
xsort=sort(x);
g=floor((percent/100)*n);

tx=xsort((g+1):(n-g));

return
