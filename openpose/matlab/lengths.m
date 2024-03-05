close all; clear all; clc;

% Reading in data
T = readtable('activity4_filtered.csv');
column_headers = {'W-TCMC','TCMC-TMCP','TMCP-TIP','TIP-TT','W-IMCP','IMCP-IPIP','IPIP-IDIP','IDIP-IT',...
    'W-MMCP','MMCP-MPIP','MPIP-MDIP','MDIP-MT','W-RMCP','RMCP-RPIP','RPIP-RDIP','RDIP-RT',...
    'W-LMCP','LMCP-LPIP','LPIP-LDIP','LDIP-LT','E-W','S-E','C-S'};
expected_length = [30,42,30,29,95,33,21,26,90,39,21,26,88,30,23,26,77,30,18,23,235,320,186];
T((end-5):end,:) = [];
time = 1:height(T);
freq = 40;

%looprange = (1:(width(angles)/3)-1)*3+1;
looprange = [4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,70,73,76];
output = [];
count = 1;
len_mean = [];
for i = looprange
    if (count <= 20) && mod(i,4) == 0
        len = sqrt((T{:,i}-T{:,76}).^2 + (T{:,i+1}-T{:,77}).^2 + (T{:,i+2}-T{:,78}).^2);
    else
        len = sqrt((T{:,i}-T{:,i-3}).^2 + (T{:,i+1}-T{:,i-2}).^2 + (T{:,i+2}-T{:,i-1}).^2);
    end
    output = cat(2,output,len);
    len_mean = [len_mean;mean(len)];
    figure
    hold on
    plot(time/freq,len)
    titlestr = strcat(column_headers{count}," length over time");
    title(titlestr)
    ylabel('length (mm)')
    xlabel('time (s)')
    set(gcf,'position',[100,300,1200,400])
    count = count + 1;
end
outtable = array2table(output);
outtable.Properties.VariableNames = column_headers;
len_err = abs(len_mean-expected_length')./(expected_length');
out_t = table(column_headers', expected_length', len_mean, len_err, 'VariableNames',["limb", "expected length (mm)", "average length (mm)", "error (%)"]);
flag = 1;