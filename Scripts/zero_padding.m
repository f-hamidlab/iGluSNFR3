%% zero_padding
% Pads data with reflected/flipped boundary values
%
% DESCRIPTION:
%   Adds padded rows at the beginning and end of a 2D array using reflected
%   boundary conditions (flips the first and last zeropad_len rows). Used to
%   handle edge cases in sliding window and filtering operations.
%
% USAGE:
%   1) data = zero_padding(data, zeropad_len)
%
% INPUTS:
%   - data: (numeric) 2D array (time x pixels/channels)
%   - zeropad_len: (int) number of rows to reflect/flip at each boundary
%
% OUTPUTS:
%   - data: (numeric) padded 2D array with reflected rows at top and bottom
%
% EXAMPLE:
%   Original: [1; 2; 3; 4; 5]
%   With zeropad_len=2, becomes: [3; 2; 1; 2; 3; 4; 5; 4; 3]
%
% Last updated: 2026-02-03 15:30

function data = zero_padding(data, zeropad_len)
    Z1 = flipud(data(1:zeropad_len, :));
    Z2 = flipud(data(end-zeropad_len+1:end,:));
    data = [Z1; data; Z2];
end
