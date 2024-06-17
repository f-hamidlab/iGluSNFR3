function peak_inds = findpeaks_2d(X, Z, MinPeakHeight)
    % Find dimensions to set up loop
    xdim = size(X,1);
    ydim = size(X,2);
    
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
end