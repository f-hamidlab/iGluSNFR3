%% zero padding
function data = zero_padding(data, zeropad_len)
    Z1 = flipud(data(1:zeropad_len, :));
    Z2 = flipud(data(end-zeropad_len+1:end,:));
    data = [Z1; data; Z2];
end
