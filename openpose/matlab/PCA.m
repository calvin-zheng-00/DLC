close all; clear all; clc;

% angles = [];
% Files = dir("C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\p2_angles");
% for i=3:length(Files)
%     filename = convertStringsToChars(Files(i).name);
%     if length(filename) < 4
%         continue
%     end
%     file = "";
%     if filename(end-3:end) == ".csv"
%         file = strcat(Files(i).folder,"\",Files(i).name);
%     else
%         continue
%     end
%     T = readtable(file);
%     T = table2array(T);
%     T = T(:,64:end);
%     angles = [angles;T];
% end

T = readtable('activity1_angles_final.csv');
angles = table2array(T);
%angles = T(:,64:end);

angles = normalize(angles,1);
[coeff,score,latent,tsquared,explained,mu] = pca(angles, 'Rows', 'pairwise');
%[coeff,score,latent,tsquared,explained,mu] = pca(angles, 'Rows', 'complete');
% coeff are the principal component coefficients, also known as loadings
%           (eigenvectors)
% score are the principal component scores
% latent are the principal component variances (eigenvalues)
% explained are the percentage of the total variance explained by each principal component
% mu are the estimated mean of each variable in X
% Name-Value arguments:
% Example: 'Algorithm','eig','Centered',false,'Rows','all','NumComponents',3
% specifies that pca uses eigenvalue decomposition algorithm, not center the data,
% use all of the observations, and return only the first three principal components.

% Find the number of components required to explain at least 95%
% variability. (idx) (set cut off two 95% and 85%)
idx = find(cumsum(explained)>95,1);
significant_85 = find(cumsum(explained)>85,1);
idx2 = find(latent>1);
var_ex = sum(explained(idx2));

% determining if two different task sets share PCAs?
% newData is the new set of data I want to apply the existing PCA onto.
%newDataPCA = (newData-mu)*coeff(:,1:idx);
% Remove idx if I want all the components

bar(explained);
title('phone typing smooth PCs percentage variance')
ylabel('Percentage explained')
xlabel('Principal Component')
writematrix(coeff, 'phone_typing_filtered_2_coeff.csv')
saveas(gcf,'Phone_typing_filtered_2.png')

% Write (specifically explained) to outout
% PC cutoff
% How to determine if two different task sets share PCAs?
% Convert pixel locations into joint angles first, or is there a way to
% combine the x,y,z coords of every point?