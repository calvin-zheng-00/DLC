close all; clear all; clc;

opt = [1,2,3,4,5;1,2,4,5,13;3,4,23,24,25;1,2,3,5,9;1,3,4,10,26;1,2,4,18,26;...
    1,2,3,4,5;1,2,3,4,10;1,3,4,6,19;2,5,18,21,26;1,2,4,5,24;1,2,3,4,5;1,2,4,5,6];
allpcarr = load("All\pconcat_coeff.mat");
allpcarr = allpcarr.coeff;

src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\PCA_coeffs';
maintemp = dir(fullfile(src,'*'));
mainfolder = {maintemp(~[maintemp.isdir]).name};
for i = 1:numel(mainfolder)
    % file = fullfile(src, mainfolder{i});
    % fprintf(1, 'Now reading %s\n', mainfolder{i});
    % thisTable = readtable(file);
    % thisTable(1,:) = [];
    theta = linspace(0,2*pi,29);
    figure
    % for j = 1:height(thisTable)
    for j = 1:5
        % rho = table2array(thisTable(:,opt(i,j)));
        rho = allpcarr(:,opt(i,j));
        rho = [rho;rho(1)];
        polarplot(theta,rho);
        hold on
    end
    thetaticks(theta*180/pi);
    thetaticklabels({'TCMC FE','TMCP FE','TIP FE','IMCP FE','IPIP FE',...
        'IDIP FE','MMCP FE','MPIP FE','MDIP FE','RMCP FE','RPIP FE',...
        'RDIP FE','LMCP FE','LPIP FE','LDIP FE','TMCP AA','TCMC R',...
        'IMCP AA','MMCP AA','RMCP AA','LMCP AA','W AA','W FE','W R',...
    	'RE FE','RS FE','RS AA','RS R'});
    rlim([-1,1]);
    set(gcf, 'Position', [100,100,700,700]);
    legend({'1','2','3','4','5'},'Location','northeastoutside','FontSize',14);
    % titlestr = strcat("PCs of grasp type: ", mainfolder{i}(1:end-10));
    % title(titlestr, 'Interpreter', 'none')
    % legend({'1','2','3','4','5','6','7','8','9','10','11','12','13','14',...
    %     '15','16','17','18','19','20','21','22','23','24','25','26','27','28'})
    hold off
end