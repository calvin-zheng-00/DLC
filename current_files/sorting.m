close all; clear all; clc;

%% Sorting trials into activities
% src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter_participant\p27_filter';
% dest = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter_participant\p27';
% 
% df = fullfile("C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/activity_order/csv/instructions_p27.csv");
% df = readtable(df,'ReadVariableNames', false);
% instruct_list = strcat("Instructions_",table2array(df(:,2)));
% instruct_no = 1;
% count = 0;
% 
% trialtemp = dir(fullfile(src,'*.csv'));
% trialfolder = {trialtemp(~[trialtemp.isdir]).name};
% for trialfolder_j = 1:numel(trialfolder)
%     data = fullfile(src,trialfolder{trialfolder_j});
%     if count > 4
%         instruct_no = instruct_no + 1;
%         count = 0;
%     end
%     outfolder = fullfile(dest,instruct_list{instruct_no});
%     if not(isfolder(outfolder))
%         mkdir(outfolder)
%     end
%     copyfile(data,outfolder)
%     count = count + 1;
% end


%% Sorting activities into categories
% src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\angles_participant';
% dest = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized_angles_2';
% 
% df = fullfile("C:\Users\czhe0008\Documents\MATLAB\Openpose\grasp_type_2.csv");
% df = readtable(df);
% grasp = dictionary(string(df.Var1),string(df.Var2));
% 
% partemp = dir(fullfile(src,'*'));
% parfolder = setdiff({partemp([partemp.isdir]).name},{'.','..'});
% for i = 1:numel(parfolder)
%     acttemp = dir(fullfile(src,parfolder{i},'*'));
%     actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
%     for j = 1:numel(actfolder)
%         % trialtemp = dir(fullfile(src,parfolder{i},actfolder{j},'*'));
%         % trialfolder = setdiff({trialtemp([trialtemp.isdir]).name},{'.','..'});
%         % for k = 1:numel(trialfolder)
%         %     file = fullfile(src,parfolder{i},actfolder{j},trialfolder{k},"df_3d.csv");
%         trialtemp = dir(fullfile(src,parfolder{i},actfolder{j},'*.csv'));
%         trialfolder = {trialtemp(~[trialtemp.isdir]).name};
%         for k = 1:numel(trialfolder)
%             file = fullfile(src,parfolder{i},actfolder{j},trialfolder{k});
%             outloc = fullfile(dest,grasp(string(actfolder{j})),trialfolder{k});
%             if not(isfolder(fullfile(dest,grasp(string(actfolder{j})))))
%                 mkdir(fullfile(dest,grasp(string(actfolder{j}))))
%             end
%             copyfile(file,outloc)
%         end
%     end
% end

%% Sorting activities together
src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter_participant';
dest = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized';

partemp = dir(fullfile(src,'*'));
parfolder = setdiff({partemp([partemp.isdir]).name},{'.','..'});
for i = 1:numel(parfolder)
    acttemp = dir(fullfile(src,parfolder{i},'*'));
    actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
    for j = 1:numel(actfolder)
        % trialtemp = dir(fullfile(src,parfolder{i},actfolder{j},'*'));
        % trialfolder = setdiff({trialtemp([trialtemp.isdir]).name},{'.','..'});
        % for k = 1:numel(trialfolder)
        %     file = fullfile(src,parfolder{i},actfolder{j},trialfolder{k},"df_3d.csv");
        trialtemp = dir(fullfile(src,parfolder{i},actfolder{j},'*.csv'));
        trialfolder = {trialtemp(~[trialtemp.isdir]).name};
        for k = 1:numel(trialfolder)
            file = fullfile(src,parfolder{i},actfolder{j},trialfolder{k});
            outloc = fullfile(dest,actfolder{j},trialfolder{k});
            if not(isfolder(fullfile(dest,actfolder{j})))
                mkdir(fullfile(dest,actfolder{j}))
            end
            copyfile(file,outloc)
        end
    end
end

%% Removing unecessary joints
% src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter_participant';
% dest = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter';
% 
% df = fullfile("C:\Users\czhe0008\Documents\MATLAB\Openpose\grasp_type_2.csv");
% df = readtable(df);
% grasp = dictionary(string(df.Var1),string(df.Var2));
% 
% partemp = dir(fullfile(src,'*'));
% parfolder = setdiff({partemp([partemp.isdir]).name},{'.','..'});
% for i = 1:numel(parfolder)
%     acttemp = dir(fullfile(src,parfolder{i},'*'));
%     actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
%     for j = 1:numel(actfolder)
%         % trialtemp = dir(fullfile(src,parfolder{i},actfolder{j},'*'));
%         % trialfolder = setdiff({trialtemp([trialtemp.isdir]).name},{'.','..'});
%         % for k = 1:numel(trialfolder)
%         %     file = fullfile(src,parfolder{i},actfolder{j},trialfolder{k},"df_3d.csv");
%         trialtemp = dir(fullfile(src,parfolder{i},actfolder{j},'*.csv'));
%         trialfolder = {trialtemp(~[trialtemp.isdir]).name};
%         for k = 1:numel(trialfolder)
%             file = fullfile(src,parfolder{i},actfolder{j},trialfolder{k});
%             df = readtable(file);
%             df(:,[1:3,79:end]) = [];
%             if not(isfolder(fullfile(dest,parfolder{i},actfolder{j})))
%                 mkdir(fullfile(dest,parfolder{i},actfolder{j}))
%             end
%             writetable(df,fullfile(dest,parfolder{i},actfolder{j},trialfolder{k}))
%         end
%     end
% end