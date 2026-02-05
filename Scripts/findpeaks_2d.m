%% findpeaks_2d
% Identifies 2D peak locations in a matrix using row and column-wise analysis
%
% DESCRIPTION:
%   Detects peaks in both x and y dimensions independently and returns indices
%   of locations that are peaks in both dimensions.
%
% USAGE:
%   1) peak_inds = findpeaks_2d(Z, MinPeakHeight)
%
% INPUTS:
%   - Z: (numeric) 2D matrix to analyze
%   - MinPeakHeight: (numeric) minimum peak height threshold
%
% OUTPUTS:
%   - peak_inds: (logical) binary matrix indicating peak locations
%
% Last updated: 2026-02-03 15:30

function peak_inds = findpeaks_2d(Z, MinPeakHeight)
    % Find dimensions to set up loop
    xdim = size(Z,1);
    ydim = size(Z,2);
    
    % Loop through x dimension to find peaks of each row
    xpeaks = zeros(size(Z));
    xwidths = NaN(size(Z));
    for i = 1:xdim
        [~,locs,w] = findpeaks(Z(i,:),'MinPeakHeight', MinPeakHeight);
        xpeaks(i,locs) = 1;
        xwidths(i,locs) = w;
    end
    
    % Loop through y dimension to find peaks of each row
    ypeaks = zeros(size(Z));
    ywidths = NaN(size(Z));
    for i = 1:ydim
        [~,locs,w] = findpeaks(Z(:,i),'MinPeakHeight', MinPeakHeight);
        ypeaks(locs,i) = 1;
        ywidths(locs,i) = w;
    end
    
    % Find indices that were peaks in both x and y
    peak_inds = xpeaks.*(xwidths>1)+ypeaks.*(ywidths>1) == 2;

    if sum(peak_inds,"all") == 0
        peak_inds = xpeaks + ypeaks == 2; % no min. requirement on width of peak
    end
end