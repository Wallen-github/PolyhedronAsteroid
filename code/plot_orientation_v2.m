%--------------------------------------------------------------------------
% Name: plot_orientation_v2
%
% Desc: Given sphericity, present a new shape model
%
% Author: Hai-Shuo Wang
% Affiliation: Univercity of Colorado Boulder, CSML
% Time: 09/26/2023
% Version 1.0:
%--------------------------------------------------------------------------
clc;
clear;
close all;
format LONG;
addpath(genpath("./data/"))

% v1: q = 1.1, v_inf = 11.1E3
info_NB = load('./data/ds_orientation_general_2b_v1.txt');
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
% beta = fliplr(beta);

figure
imagesc(alpha,beta,info_matrix);
colormap("pink");
colorbar;
title('relative shift distance (m)');
xlabel('\alpha (deg)');
ylabel('\beta (deg)');
xticks([0 20 40 60 80 100 120 140 160 180]);
xticklabels({'0','20','40','60','80','100','120','140','160','180'});
yticks([0 40 80 120 160 200 240 280 320 360]);
yticklabels({'360','320','280','240','200','160','120','80','40','0'});
% colorbar(ticks=[0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6],TickLabels={'disrupt','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6'});
set(gca,Fontsize=20)


figure
imagesc(alpha,beta,Pe_ratio);
colormap("pink");
colorbar;
title('Period Ratio (P_l/P_0)');
xlabel('\alpha (deg)');
ylabel('\beta (deg)');
xticks([0 20 40 60 80 100 120 140 160 180]);
xticklabels({'0','20','40','60','80','100','120','140','160','180'});
yticks([0 40 80 120 160 200 240 280 320 360]);
yticklabels({'360','320','280','240','200','160','120','80','40','0'});
% ds_orbit_4b_v1:
% colorbar(ticks=[0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6],TickLabels={'disrupt','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6'});
set(gca,Fontsize=20)