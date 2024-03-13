% takes an array of data (assuming each row is a different channel)
% and resamples each row
% returns the resampled array with rows in the same order
%
% inputs: 
% x - data for resampling
% fsnew - new sampling rate
% fs - original sampling rate
%
% outputs: 
% y - resampled data
% Credit to Hannah Tay

function y = resample_array(x, fsnew, fs)
y = [];
    for i = 1:size(x,1)
        row = x(i,:);
        resampled_row = resample(row, fsnew, fs);
        y = [y; resampled_row];
    end
end