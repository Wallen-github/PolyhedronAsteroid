%--------------------------------------------------------------------------
% Name: plot_pervelNB_v3
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

% Ig_orbit_2b_v1: vel=[0,20]E3, perigee=[1.1,8.3]RE, surface 0.005
% Ig_orbit_2b_v2: vel=[0,20]E3, perigee=[1.1,8.3]RE, surface 0.0005
% Ig_orbit_2b_v3: vel=[0,20]E3, perigee=[1.1,2.5]RE, surface 0.005
% Ig_orbit_2b_v4: vel=[0,20]E3, perigee=[1.1,5.0]RE, surface 0.005, Ig-I0
% Ig_orbit_2b_v5: vel=[0,20]E3, perigee=[1.1,5.0]RE, surface 0.00005, Ig-I0
% Ig_orbit_3b_v1: vel=[0,20]E3, perigee=[1.1,8.3]RE
% Ig_orbit_3b_v2: vel=[0,20]E3, perigee=[1.1,2.5]RE
% Ig_orbit_3b_v3: vel=[0,20]E3, perigee=[1.1,5.0]RE, Ig-I0
% Ig_orbit_4b_v1: vel=[0,20]E3, perigee=[1.1,8.3]RE
% Ig_orbit_4b_v2: vel=[0,20]E3, perigee=[1.1,2.5]RE
% Ig_orbit_4b_v3: vel=[0,20]E3, perigee=[1.1,5.0]RE, Ig-I0

% ds_orbit_2b_v1: vel=[0,20]E3, perigee=[1.1,6.0]RE, surface 0.005
% ds_orbit_2b_v2: vel=[0,20]E3, perigee=[1.1,6.0]RE, surface 0.00005
% ds_orbit_3b_v1: vel=[0,20]E3, perigee=[1.1,6.0]RE
% ds_orbit_4b_v1: vel=[0,20]E3, perigee=[1.1,6.0]RE
% ds_orientation_Apophis_2b_v1: alpha=[0,180], beta=[0,360], surface 0.005 
info_NB = load('./data/ds_orbit_4b_v1.txt');
nbR3 = 2;
row_num = 10;
col_num = 10;

for i=1:col_num
    
    per(i) = info_NB(i + (i-1)*row_num,1);
    vel(i) = info_NB(i,2);

    for j=1:row_num

        Id_ratio(j,i) = info_NB(j + (i-1)*col_num,3);
        Pe_ratio(j,i) = info_NB(j + (i-1)*col_num,4);
        Mo_ratio(j,i) = info_NB(j + (i-1)*col_num,5);
        Bre_info(j,i) = info_NB(j + (i-1)*col_num,6);
        
        info_matrix(j,i) = log10(Id_ratio(j,i)*168);
%         info_matrix(j,i) = Id_ratio(j,i)*168;
        if Bre_info(j,i) == 0
            info_matrix(j,i) = nan;
            Pe_ratio(j,i) = nan;
        end

    end
end

Id_ratio = flipud(Id_ratio);
% Pe_ratio = flipud(Pe_ratio);
Mo_ratio = flipud(Mo_ratio);
% info_matrix = flipud(info_matrix);
vel = flipud(vel');

figure
imagesc(per,vel,info_matrix);
colormap("pink");
colorbar;
title('relative shift distance (meter)');
xlabel('perigee (R_E)');
ylabel('velocity_{\infty} (km/s)');
xticks([1.1 1.6 2.2 2.7 3.3 3.8 4.4 4.9 5.5 6]);
xticklabels({'1.1','1.6','2.2','2.7','3.3','3.8','4.4','4.9','5.5','6'});
yticks([0 2.2 4.4 6.7 8.9 11.1 13.3 15.6 17.8 20]);
yticklabels({'20','17.8','15.6','13.3','11.1','8.9','6.7','4.4','2.2','0'});
% ds_orbit_3b_v1: ds_orbit_4b_v1:
colorbar(ticks=[-1.7,-1.5,-1,-0.5,0,0.5,1],TickLabels={'disrupt','0.03','0.1','0.3','1.0','3.1','10'});
% ds_orbit_2b_v2:
% colorbar(ticks=[0.8,0.81,0.82,0.83,0.84,0.85,0.86],TickLabels={'disrupt','6.4','6.6','6.7','6.9','7.0','7.2'});
% ds_orbit_2b_v1:
% colorbar(ticks=[-4.9,-4,-3,-2,-1,0],TickLabels={'disrupt','10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'});
set(gca,Fontsize=20)

figure
imagesc(per,vel,Pe_ratio);
colormap("pink");
colorbar;
title('Period Ratio (P_l/P_0)');
xlabel('perigee (R_E)');
ylabel('velocity_{\infty} (km/s)');
xticks([1.1 1.6 2.2 2.7 3.3 3.8 4.4 4.9 5.5 6]);
xticklabels({'1.1','1.6','2.2','2.7','3.3','3.8','4.4','4.9','5.5','6'});
yticks([0 2.2 4.4 6.7 8.9 11.1 13.3 15.6 17.8 20]);
yticklabels({'20','17.8','15.6','13.3','11.1','8.9','6.7','4.4','2.2','0'});
% ds_orbit_4b_v1:
colorbar(ticks=[0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6],TickLabels={'disrupt','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6'});
% ds_orbit_3b_v1:
% colorbar(ticks=[0.1, 0.2, 0.6, 1.0, 1.4],TickLabels={'disrupt','0.2','0.6','1.0','1.4'});
% ds_orbit_2b_v2:
% colorbar(ticks=[0.5, 1, 3, 5, 7, 9],TickLabels={'disrupt','1','3','5','7','9'});
% ds_orbit_2b_v1:
% colorbar(ticks=[1,10,20,30,40],TickLabels={'disrupt','10','20','30','40'});
set(gca,Fontsize=20)

