%% remv_zero_padding
% Removes reflected padding added by zero_padding function
%
% DESCRIPTION:
%   Removes the padded rows from the beginning and end of a 2D array,
%   reversing the effect of zero_padding. Extracts the central portion of
%   a padded array to return to original dimensions.
%
% USAGE:
%   1) data = remv_zero_padding(data, zeropad_len)
%
% INPUTS:
%   - data: (numeric) padded 2D array (time x pixels/channels)
%   - zeropad_len: (int) number of rows that were padded at each boundary
%
% OUTPUTS:
%   - data: (numeric) unpadded 2D array with central portion extracted
%
% EXAMPLE:
%   Padded: [3; 2; 1; 2; 3; 4; 5; 4; 3]
%   With zeropad_len=2, returns: [1; 2; 3; 4; 5]
%
% Last updated: 2026-02-03 15:30

function data = remv_zero_padding(data, zeropad_len)
    
    data = data(zeropad_len+1:size(data,1)-zeropad_len,:);
end