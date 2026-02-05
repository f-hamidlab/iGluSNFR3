%% knee_pt
% Returns the x-location of a knee/elbow point of a curve
%
% DESCRIPTION:
%   Finds the knee point of a curve y=f(x) by fitting two lines (one to the left
%   and one to the right of each bisection point) and finding the point that
%   minimizes the sum of fitting errors. Useful for determining cutoff values
%   like eigenvalue decay or elbow method in clustering.
%
% USAGE:
%   1) [res_x, idx_of_result] = knee_pt(y)
%   2) [res_x, idx_of_result] = knee_pt(y, x)
%   3) [res_x, idx_of_result] = knee_pt(y, x, just_return)
%
% INPUTS:
%   - y: (numeric) vector of data values, length >= 3 elements (required)
%   - x: (numeric) x-coordinates corresponding to y. If not provided, defaults to 1:length(y)
%   - just_return: (logical) if true, returns NaN on error instead of erroring out
%
% OUTPUTS:
%   - res_x: (numeric) x-coordinate at the knee point
%   - idx_of_result: (numeric) index of the knee point in the x vector
%
% NOTES:
%   - Function will never return the first or last point as the knee
%   - x and y do not need to be sorted, but must correspond
%   - Errors are calculated as sum(abs(residuals)) by default
%
% Last updated: 2026-02-03 15:30

function [res_x, idx_of_result] = knee_pt(y,x,just_return)
% 
%dk /2012
%{
% test vectors:
[i,ix]=knee_pt([30:-3:12,10:-2:0])  %should be 7 and 7
knee_pt([30:-3:12,10:-2:0]')  %should be 7
knee_pt(rand(3,3))  %should error out
knee_pt(rand(3,3),[],false)  %should error out
knee_pt(rand(3,3),[],true)  %should return Nan
knee_pt([30:-3:12,10:-2:0],[1:13])  %should be 7
knee_pt([30:-3:12,10:-2:0],[1:13]*20)  %should be 140
knee_pt([30:-3:12,10:-2:0]+rand(1,13)/10,[1:13]*20)  %should be 140
knee_pt([30:-3:12,10:-2:0]+rand(1,13)/10,[1:13]*20+rand(1,13)) %should be close to 140
x = 0:.01:pi/2; y = sin(x); [i,ix]=knee_pt(y,x)  %should be around .9 andaround 90
[~,reorder]=sort(rand(size(x)));xr = x(reorder); yr=y(reorder);[i,ix]=knee_pt(yr,xr)  %i should be the same as above and xr(ix) should be .91
knee_pt([10:-1:1])  %degenerate condition -- returns location of the first "knee" error minimum: 2
%}
%set internal operation flags
use_absolute_dev_p = true;  %ow quadratic
%deal with issuing or not not issuing errors
issue_errors_p = true;
if (nargin > 2 && ~isempty(just_return) && just_return)
    issue_errors_p = false;
end
%default answers
res_x = nan;
idx_of_result = nan;
%check...
if (isempty(y))
    if (issue_errors_p)
        error('knee_pt: y can not be an empty vector');
    end
    return;
end
%another check
if (sum(size(y)==1)~=1)
    if (issue_errors_p)
        error('knee_pt: y must be a vector');
    end
    
    return;
end
%make a vector
y = y(:);
%make or read x
if (nargin < 2 || isempty(x))
    x = (1:length(y))';
else
    x = x(:);
end
%more checking
if (ndims(x)~= ndims(y) || ~all(size(x) == size(y)))
    if (issue_errors_p)
        error('knee_pt: y and x must have the same dimensions');
    end
    
    return;
end
%and more checking
if (length(y) < 3)
    if (issue_errors_p)
        error('knee_pt: y must be at least 3 elements long');
    end
    return;
end
%make sure the x and y are sorted in increasing X-order
if (nargin > 1 && any(diff(x)<0))
    [~,idx]=sort(x);
    y = y(idx);
    x = x(idx);
else
    idx = 1:length(x);
end
%the code below "unwraps" the repeated regress(y,x) calls.  It's
%significantly faster than the former for longer y's
%
%figure out the m and b (in the y=mx+b sense) for the "left-of-knee"
sigma_xy = cumsum(x.*y);
sigma_x  = cumsum(x);
sigma_y  = cumsum(y);
sigma_xx = cumsum(x.*x);
n        = (1:length(y))';
det = n.*sigma_xx-sigma_x.*sigma_x;
mfwd = (n.*sigma_xy-sigma_x.*sigma_y)./det;
bfwd = -(sigma_x.*sigma_xy-sigma_xx.*sigma_y) ./det;
%figure out the m and b (in the y=mx+b sense) for the "right-of-knee"
sigma_xy = cumsum(x(end:-1:1).*y(end:-1:1));
sigma_x  = cumsum(x(end:-1:1));
sigma_y  = cumsum(y(end:-1:1));
sigma_xx = cumsum(x(end:-1:1).*x(end:-1:1));
n        = (1:length(y))';
det = n.*sigma_xx-sigma_x.*sigma_x;
mbck = flipud((n.*sigma_xy-sigma_x.*sigma_y)./det);
bbck = flipud(-(sigma_x.*sigma_xy-sigma_xx.*sigma_y) ./det);
%figure out the sum of per-point errors for left- and right- of-knee fits
error_curve = nan(size(y));
for breakpt = 2:length(y-1)
    delsfwd = (mfwd(breakpt).*x(1:breakpt)+bfwd(breakpt))-y(1:breakpt);
    delsbck = (mbck(breakpt).*x(breakpt:end)+bbck(breakpt))-y(breakpt:end);
    %disp([sum(abs(delsfwd))/length(delsfwd), sum(abs(delsbck))/length(delsbck)])
    if (use_absolute_dev_p)
        % error_curve(breakpt) = sum(abs(delsfwd))/sqrt(length(delsfwd)) + sum(abs(delsbck))/sqrt(length(delsbck));
        error_curve(breakpt) = sum(abs(delsfwd))+ sum(abs(delsbck));
    else
        error_curve(breakpt) = sqrt(sum(delsfwd.*delsfwd)) + sqrt(sum(delsbck.*delsbck));
    end
end
%find location of the min of the error curve
[~,loc] = min(error_curve);
res_x = x(loc);
idx_of_result = idx(loc);
end
