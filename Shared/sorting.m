close all; clear all; clc;

%% Sorting trials into activities
% src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\angles\p24';
% dest = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\angles';
% 
% df = fullfile("C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/activity_order/csv/instructions_p24.csv");
% df = readmatrix(df);
% instruct_list = df.videos;
% instruct_no = 0;
% count = -1;
% 
% trialtemp = dir(fullfile(src,'*.csv'));
% trialfolder = {trialtemp(~[trialtemp.isdir]).name};
% for trialfolder_j = 1:numel(trialfolder)
%     data = fullfile(src,catfolder{catfolder_i},trialfolder{trialfolder_j});
%     outfolder = fullfile(dest,catfolder{catfolder_i});
%     outloc = fullfile(outfolder,strcat(trialfolder{trialfolder_j}(1:end-13),"_angles.csv"));
% end


%% Sorting activities into categories
src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\3d';
dest = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized_2';

df = fullfile("C:\Users\czhe0008\Documents\MATLAB\Openpose\grasp_type_2.csv");
df = readtable(df);
grasp = dictionary(string(df.Var1),string(df.Var2));

partemp = dir(fullfile(src,'*'));
parfolder = setdiff({partemp([partemp.isdir]).name},{'.','..'});
for i = 1:numel(parfolder)
    acttemp = dir(fullfile(src,parfolder{i},'*'));
    actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
    for j = 1:numel(actfolder)
        trialtemp = dir(fullfile(src,parfolder{i},actfolder{j},'*'));
        trialfolder = setdiff({trialtemp([trialtemp.isdir]).name},{'.','..'});
        for k = 1:numel(trialfolder)
            file = fullfile(src,parfolder{i},actfolder{j},trialfolder{k},"df_3d.csv");
            if exist(file, 'file')
                outloc = fullfile(dest,grasp(string(actfolder{j})),strcat(trialfolder{k},".csv"));
                if not(isfolder(fullfile(dest,grasp(string(actfolder{j})))))
                    mkdir(fullfile(dest,grasp(string(actfolder{j}))))
                end
                copyfile(file,outloc)
            else
                fprintf('Warning: file does not exist:\n%s', file);
            end
        end
    end
end