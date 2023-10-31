function data = remv_zero_padding(data, zeropad_len)
    
    data = data(zeropad_len+1:size(data,1)-zeropad_len,:);
end