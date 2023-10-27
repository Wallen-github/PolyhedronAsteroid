%--------------------------------------------------------------------------
% Name: plot_orientation
%
% Desc: Given sphericity, present a new shape model
%
% Author: Hai-Shuo Wang
% Affiliation: Univercity of Colorado Boulder, CSML
% Time: 09/21/2023
% Version 1.0:
%--------------------------------------------------------------------------
clc;
clear;
close all;
format LONG;
addpath(genpath("./data/"))

% v1: q = 1.1, v_inf = 11.1E3
info_NB = load('./data/ds_orientation_apophis_2b_v1.txt');
nbR3 = 2;
row_num = 10;
col_num = 10;

for i=1:col_num
    
    alpha(i) = info_NB(i + (i-1)*row_num,1);
    beta(i) = info_NB(i,2);

    for j=1:row_num

        Id_ratio(j,i) = info_NB(j + (i-1)*col_num,3);
        Pe_ratio(j,i) = info_NB(j + (i-1)*col_num,4);
        Mo_ratio(j,i) = info_NB(j + (i-1)*col_num,5);
        Bre_info(j,i) = info_NB(j + (i-1)*col_num,6);
        
        info_matrix(j,i) = Id_ratio(j,i)*168;
        if Bre_info(j,i) == 0
            info_matrix(j,i) = nan;
            Pe_ratio(j,i) = nan;
        end

    end
end

Id_ratio = flipud(Id_ratio);
Pe_ratio = flipud(Pe_ratio);
Mo_ratio = flipud(Mo_ratio);
info_matrix = flipud(info_matrix);
beta = fliplr(beta);

figure
h = heatmap(alpha,beta,info_matrix,Colormap=pink,CellLabelColor='none');
h.Title = 'relative shift distance (m)';
h.XLabel = 'alpha (deg)';
h.YLabel = 'beta (deg)';
set(gca,Fontsize=20)

figure
h = heatmap(alpha,beta,Pe_ratio,Colormap=pink,CellLabelColor='none');
h.Title = 'Period Ratio (P_l/P_0)';
h.XLabel = 'alpha (deg)';
h.YLabel = 'beta (deg)';
set(gca,Fontsize=20)