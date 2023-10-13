close all; clear all; clc;

% Reading in data
File = readtable('eating_soup.csv');

T=File;
time = 7;

% Cleaning data so that only the x,y,z values remain
T(:,(end-13):end) = [];  % removing extra stats at the end of the table
T(:,4:6:end) = []; % removing error values
T(:,4:5:end) = []; % removing number of cams
T(:,4:4:end) = []; % removing score values

A = table2array(T);
[coeff,score,latent,tsquared,explained,mu] = pca(A, 'Rows', 'pairwise');
% coeff are the principal component coefficients, also known as loadings
% score are the principal component scores
% latent are the principal component variances
% explained are the percentage of the total variance explained by each principal component
% mu are the estimated mean of each variable in X

% Name-Value arguments:
% Example: 'Algorithm','eig','Centered',false,'Rows','all','NumComponents',3
% specifies that pca uses eigenvalue decomposition algorithm, not center the data,
% use all of the observations, and return only the first three principal components.


% Find the number of components required to explain at least 95% variability.
idx = find(cumsum(explained)>95,1);
% newData is the new set of data I want to apply the existing PCA onto.
newDataPCA = (newData-mu)*coeff(:,1:idx);
% Remove idx if I want all the components

% Loop through multiple files
% Write (specifically explained) to outout
% PC cutoff
% How to determine if two different task sets share PCAs?
% Convert pixel locations into joint angles first, or is there a way to
% combine the x,y,z coords of every point?